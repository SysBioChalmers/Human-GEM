"""Task files are checked with the rules RAVEN's checkTasks applies before solving."""

import cobra

import testMetabolicTasks

HEADER = "\tID\tDESCRIPTION\tSHOULD FAIL\tIN\tIN LB\tIN UB\tOUT\tOUT LB\tOUT UB\n"


def _model():
    model = cobra.Model("m")
    model.compartments = {"c": "cytosol", "e": "extracellular"}
    reaction = cobra.Reaction("R1")
    reaction.add_metabolites({
        cobra.Metabolite("h_c", name="H+", compartment="c"): -1,
        cobra.Metabolite("glc_e", name="glucose", compartment="e"): 1,
    })
    model.add_reactions([reaction])
    return model


def _problems(tmp_path, rows):
    path = tmp_path / "metabolicTasks_Test.txt"
    path.write_text(HEADER + "".join("\t" + "\t".join(row) + "\n" for row in rows), encoding="utf-8")
    return testMetabolicTasks.task_file_problems(_model(), path)


def test_valid_task_has_no_problems(tmp_path):
    rows = [["1", "ok", "", "glucose[e];ALLMETSIN[e]", "", "", "H+[c]", "", ""],
            ["", "", "", "H+[c]", "", "", "ALLMETS", "", ""]]
    assert _problems(tmp_path, rows) == []


def test_unknown_metabolite_is_reported(tmp_path):
    problems = _problems(tmp_path, [["1", "t", "", "fructose[e]", "", "", "H+[m]", "", ""]])
    assert len(problems) == 2
    assert "'fructose[e]' is not in the model" in problems[0]
    assert "'H+[m]' is not in the model" in problems[1]


def test_input_and_output_with_nonzero_lower_bounds_is_reported(tmp_path):
    # CellfieConsensus task 21: H+[c] IN 1-1 and OUT 1-1.
    problems = _problems(tmp_path, [["21", "t", "", "H+[c]", "1", "1", "H+[c]", "1", "1"]])
    expected = "metabolicTasks_Test.txt [21] t: H+[C] is an input and an output, both with a nonzero lower bound"
    assert problems == [expected]


def test_input_and_output_with_one_zero_lower_bound_is_accepted(tmp_path):
    assert _problems(tmp_path, [["1", "t", "", "H+[c]", "", "", "H+[c]", "1", "1"]]) == []


def test_compartment_case_is_ignored_as_in_raven(tmp_path):
    rows = [["1", "upper-case compartment", "", "glucose[E]", "", "", "H+[C]", "", ""]]
    assert _problems(tmp_path, rows) == []
