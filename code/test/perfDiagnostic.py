"""Isolated performance diagnostic for the gene-essentiality CI slowdown.

On 2026-09-10 (run 34494319341, raven-toolbox @ 61da432), building DLD1's context
model took 23 minutes. Re-running that *exact* code and raven-toolbox commit the
next day (run 34645940348 attempt 1, same 61da432) took over 5 hours and did not
finish. No code changed between those two runs, which rules out anything in
Human-GEM or raven-toolbox as the cause -- something about the execution
environment (Gurobi WLS server throughput, or the GitHub-hosted runner's actual
delivered CPU) got much worse in between.

This script isolates each candidate cause with a fast (~5 min), self-contained
benchmark that has no dependency on cobra, raven-toolbox, or the Human-GEM model:

  * CPU: a fixed-size numpy matmul, unrelated to Gurobi or any license, to
    measure raw floating-point throughput independent of solver/licensing.
  * Gurobi/WLS handshake: time to acquire a licensed environment, separate
    from any solve.
  * Gurobi MILP throughput: solve fixed-seed random 0/1 multi-knapsack MIPs of
    increasing size, all well under any time limit, so runtime measures
    solver throughput and not "did it time out." Gurobi's own reported
    NodeCount/IterCount/Runtime give a throughput figure that isn't confounded
    by wall-clock scheduling noise the way an external timer alone would be.

None of this reproduces ftINIT's specific problem structure -- it is a
proxy for "is this environment's raw compute/licensing throughput degraded",
not a replacement for the real workload.
"""

from __future__ import annotations

import platform
import time

import numpy as np


def _log(msg: str) -> None:
    print(f"[perf-diag] {msg}", flush=True)


def cpu_info() -> None:
    _log(f"platform: {platform.platform()}")
    _log(f"processor: {platform.processor()!r}")
    try:
        with open("/proc/cpuinfo") as fh:
            text = fh.read()
        model = next((l for l in text.splitlines() if l.startswith("model name")), None)
        mhz = [l for l in text.splitlines() if l.startswith("cpu MHz")]
        count = text.count("processor\t:")
        _log(f"/proc/cpuinfo: {count} logical CPU(s); {model or 'model name unavailable'}")
        if mhz:
            _log(f"/proc/cpuinfo cpu MHz (first 4): {mhz[:4]}")
    except OSError as exc:
        _log(f"/proc/cpuinfo unavailable: {exc}")
    import os

    _log(f"os.cpu_count(): {os.cpu_count()}")
    try:
        quota_path = "/sys/fs/cgroup/cpu.max"
        with open(quota_path) as fh:
            _log(f"cgroup cpu.max: {fh.read().strip()}")
    except OSError:
        pass


def cpu_benchmark(n: int = 2000, repeats: int = 3) -> None:
    _log(f"--- CPU benchmark: {n}x{n} float64 matmul, {repeats} repeat(s) ---")
    rng = np.random.default_rng(0)
    a = rng.standard_normal((n, n))
    b = rng.standard_normal((n, n))
    times = []
    for i in range(repeats):
        t0 = time.perf_counter()
        c = a @ b
        elapsed = time.perf_counter() - t0
        times.append(elapsed)
        gflops = 2 * n**3 / elapsed / 1e9
        _log(f"  run {i + 1}/{repeats}: {elapsed:.3f}s, {gflops:.2f} GFLOP/s (checksum {c[0, 0]:.4f})")
    _log(f"  median: {sorted(times)[len(times) // 2]:.3f}s")


def gurobi_handshake() -> "gurobipy.Env | None":
    import gurobipy

    _log("--- Gurobi WLS handshake ---")
    t0 = time.perf_counter()
    try:
        env = gurobipy.Env()
    except gurobipy.GurobiError as exc:
        _log(f"  FAILED to acquire environment: {exc}")
        return None
    elapsed = time.perf_counter() - t0
    _log(f"  environment acquired in {elapsed:.2f}s")
    return env


def _random_knapsack(rng: np.random.Generator, n_vars: int, n_constrs: int):
    weights = rng.integers(1, 100, size=(n_constrs, n_vars))
    capacities = weights.sum(axis=1) * 0.5
    values = rng.integers(1, 100, size=n_vars)
    return weights, capacities, values


def gurobi_mip_benchmark(env, sizes=((80, 8), (150, 12), (220, 15)), time_limit=120.0) -> None:
    import gurobipy
    from gurobipy import GRB

    _log(f"--- Gurobi MILP throughput: {len(sizes)} fixed-seed multi-knapsack instance(s) ---")
    for idx, (n_vars, n_constrs) in enumerate(sizes):
        rng = np.random.default_rng(1000 + idx)
        weights, capacities, values = _random_knapsack(rng, n_vars, n_constrs)

        model = gurobipy.Model(env=env)
        model.Params.OutputFlag = 0
        model.Params.TimeLimit = time_limit
        model.Params.Seed = 0
        x = model.addVars(n_vars, vtype=GRB.BINARY, name="x")
        for c in range(n_constrs):
            model.addConstr(gurobipy.quicksum(weights[c, j] * x[j] for j in range(n_vars)) <= capacities[c])
        model.setObjective(gurobipy.quicksum(values[j] * x[j] for j in range(n_vars)), GRB.MAXIMIZE)

        t0 = time.perf_counter()
        model.optimize()
        wall = time.perf_counter() - t0

        status = model.Status
        nodes = model.NodeCount
        iters = model.IterCount
        runtime = model.Runtime
        gap = model.MIPGap if status == GRB.OPTIMAL else float("nan")
        node_rate = nodes / runtime if runtime > 0 else float("nan")
        iter_rate = iters / runtime if runtime > 0 else float("nan")
        _log(
            f"  [{n_vars}vars/{n_constrs}constrs] status={status} wall={wall:.2f}s "
            f"gurobi_runtime={runtime:.2f}s nodes={nodes:.0f} ({node_rate:.1f}/s) "
            f"simplex_iters={iters:.0f} ({iter_rate:.1f}/s) gap={gap}"
        )
        if status == GRB.TIME_LIMIT:
            _log("  WARNING: hit time_limit without proving optimality -- instance too hard for this run, "
                 "throughput figures above are a floor, not a clean rate")


def main() -> None:
    cpu_info()
    cpu_benchmark()
    env = gurobi_handshake()
    if env is not None:
        gurobi_mip_benchmark(env)
        env.dispose()
    else:
        _log("Skipping MILP benchmark: no Gurobi environment")


if __name__ == "__main__":
    main()
