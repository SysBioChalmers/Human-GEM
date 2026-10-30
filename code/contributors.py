"""Keep the Human-GEM contributor list: `.all-contributorsrc` and the README table.

Replaces the all-contributors CLI, so that crediting needs Python only and can be
run and tested locally. The README table is written in the same format as the CLI
writes it.

Usage (from the repository root):
    python code/contributors.py add <login> <type>[,<type>...]
    python code/contributors.py generate
    python code/contributors.py check

`add` reads the name, avatar and profile link from the GitHub API, using the token
in GH_TOKEN or GITHUB_TOKEN if one is set. The commands `apply`, `suggest` and
`thank` are run by `.github/workflows/add-contributor.yml`.
"""

import argparse
import json
import os
import re
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RC = ROOT / ".all-contributorsrc"
README = ROOT / "README.md"
REPO = "SysBioChalmers/Human-GEM"

# The six types Human-GEM uses, in the order they are shown. Each maps to its
# emoji, its title and the link the all-contributors CLI gives it ("{login}" is
# filled in; None links to an anchor in the README).
TYPES = {
    "bug": ("🐛", "Bug reports", f"https://github.com/{REPO}/issues?q=author%3A{{login}}"),
    "research": ("🔬", "Research", None),
    "code": ("💻", "Code", f"https://github.com/{REPO}/commits?author={{login}}"),
    "review": ("👀", "Reviewed Pull Requests",
               f"https://github.com/{REPO}/pulls?q=is%3Apr+reviewed-by%3A{{login}}"),
    "doc": ("📖", "Documentation", f"https://github.com/{REPO}/commits?author={{login}}"),
    "ideas": ("🤔", "Ideas, Planning, & Feedback", None),
}
# Words people use for a type, mapped to the type.
ALIASES = {"bugs": "bug", "docs": "doc", "documentation": "doc", "idea": "ideas",
           "reviews": "review", "data": "research", "evidence": "research"}

LOGIN = r"[A-Za-z0-9](?:[A-Za-z0-9-]{0,38})"
# "/credit @user bug, code", one per line. The older
# "@all-contributors please add @user for bug, code" is still understood.
COMMAND = re.compile(
    rf"^\s*(?:/credit\s+@?(?P<u1>{LOGIN})(?:\s+(?:for\s+)?(?P<t1>.*?))?"
    rf"|@all-contributors\s+(?:please\s+)?add\s+@?(?P<u2>{LOGIN})\s+for\s+(?P<t2>.*?))"
    r"[.!]?\s*$",
    re.IGNORECASE | re.MULTILINE,
)

LIST_START = "<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->"
LIST_END = "<!-- ALL-CONTRIBUTORS-LIST:END -->"
BADGE = re.compile(r"(<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->\n)"
                   r".*?(\n<!-- ALL-CONTRIBUTORS-BADGE:END -->)", re.DOTALL)


def load():
    return json.loads(RC.read_text(encoding="utf-8"))


def save(rc):
    text = json.dumps(rc, indent=2, ensure_ascii=False) + "\n"
    RC.write_bytes(text.replace("\n", newline_of(RC)).encode("utf-8"))


def newline_of(path):
    return "\r\n" if b"\r\n" in path.read_bytes() else "\n"


def find(rc, login):
    return next((c for c in rc["contributors"] if c["login"].lower() == login.lower()), None)


def parse_types(text):
    """Split "bug, code and doc" into known types; return (types, unknown words)."""
    words = [w.strip().lower() for w in re.split(r",|\s+and\s+|\s+", text) if w.strip()]
    types, unknown = [], []
    for w in words:
        t = ALIASES.get(w, w)
        if t in TYPES:
            if t not in types:
                types.append(t)
        else:
            unknown.append(w)
    return types, unknown


def parse(body):
    """Read every credit command in a comment.

    Returns (credits, errors): credits is a list of (login, types), errors a list
    of messages for commands that name no known type.
    """
    credits, errors = [], []
    for m in COMMAND.finditer(body or ""):
        login = m.group("u1") or m.group("u2")
        types, unknown = parse_types((m.group("t1") if m.group("u1") else m.group("t2")) or "")
        if unknown:
            errors.append(f"`@{login}`: unknown type {', '.join(f'`{w}`' for w in unknown)}")
        elif not types:
            errors.append(f"`@{login}`: no contribution type given")
        else:
            credits.append((login, types))
    return credits, errors


def github_user(login):
    req = urllib.request.Request(f"https://api.github.com/users/{login}",
                                 headers={"Accept": "application/vnd.github+json"})
    token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    with urllib.request.urlopen(req) as r:
        return json.load(r)


def add(rc, login, types, user=None):
    """Credit `login` for `types`; return the types that were new for them."""
    entry = find(rc, login)
    if entry is None:
        user = user or github_user(login)
        blog = (user.get("blog") or "").strip()
        if blog and not blog.startswith("http"):
            blog = "https://" + blog
        entry = {
            "login": user["login"],
            "name": user.get("name") or user["login"],
            "avatar_url": user["avatar_url"],
            "profile": blog or user["html_url"],
            "contributions": [],
        }
        rc["contributors"].append(entry)
    new = [t for t in types if t not in entry["contributions"]]
    entry["contributions"] = [t for t in TYPES if t in entry["contributions"] + new]
    return new


def cell(c):
    links = []
    for t in c["contributions"]:
        emoji, title, url = TYPES[t]
        href = url.format(login=c["login"]) if url else f"#{t}-{c['login']}"
        links.append(f'<a href="{href}" title="{title}">{emoji}</a>')
    return (f'      <td align="center" valign="top" width="10%"><a href="{c["profile"]}">'
            f'<img src="{c["avatar_url"]}?s=64" width="64px;" alt="{c["name"]}"/><br />'
            f'<sub><b>{c["name"]}</b></sub></a><br />{" ".join(links)}</td>')


def table(rc):
    people = sorted(rc["contributors"], key=lambda c: c["name"])
    per = rc.get("contributorsPerLine", 10)
    rows = []
    for i in range(0, len(people), per):
        rows += ["    <tr>", *(cell(c) for c in people[i:i + per]), "    </tr>"]
    return "\n".join([LIST_START, "<!-- prettier-ignore-start -->", "<!-- markdownlint-disable -->",
                      "<table>", "  <tbody>", *rows, "  </tbody>", "</table>",
                      "<!-- markdownlint-restore -->", "<!-- prettier-ignore-end -->", LIST_END])


def render(rc, readme):
    """Return the README text with the badge and the table rebuilt from `rc`."""
    nl = "\r\n" if "\r\n" in readme else "\n"
    text = readme.replace("\r\n", "\n")
    start, end = text.index(LIST_START), text.index(LIST_END) + len(LIST_END)
    text = text[:start] + table(rc) + text[end:]
    badge = rc["badgeTemplate"].replace("<%= contributors.length %>", str(len(rc["contributors"])))
    text = BADGE.sub(lambda m: m.group(1) + badge + m.group(2), text)
    return text.replace("\n", nl)


def generate(rc):
    README.write_bytes(render(rc, README.read_bytes().decode("utf-8")).encode("utf-8"))


def problems(rc, readme):
    """Errors in the list: unknown types, repeated logins, a README out of date."""
    out, seen = [], set()
    for c in rc["contributors"]:
        if c["login"].lower() in seen:
            out.append(f"{c['login']} is listed more than once")
        seen.add(c["login"].lower())
        bad = [t for t in c["contributions"] if t not in TYPES]
        if bad:
            out.append(f"{c['login']} has unknown types {bad}; use {', '.join(TYPES)}")
        if not c["contributions"]:
            out.append(f"{c['login']} has no contribution type")
    if not out and render(rc, readme) != readme:
        out.append("README.md does not match .all-contributorsrc; run `python code/contributors.py generate`")
    return out


def missing(rc, logins):
    """The logins that are not in the list, without bots, in their first order."""
    out = []
    for login in logins:
        if login and not login.endswith("[bot]") and not find(rc, login) and login not in out:
            out.append(login)
    return out


def suggest(rc, pr_author, pr_files, issues):
    """Proposed credits for a merged pull request, for people not in the list yet.

    `issues` is a list of (author, labels) of the issues the pull request closes.
    The author of an issue is proposed for `ideas` if it is labelled `feature`
    and `bug` otherwise. The author of the pull request is proposed for `research`
    for changes to model/ or data/, `code` for code, and `doc` for Markdown.
    """
    types = {}
    for author, labels in issues:
        types.setdefault(author, set()).add("ideas" if "feature" in labels else "bug")
    pr = types.setdefault(pr_author, set())
    for f in pr_files:
        if f.startswith(("model/", "data/")):
            pr.add("research")
        elif f.endswith(".md"):
            pr.add("doc")
        else:
            pr.add("code")
    return [(login, [t for t in TYPES if t in types[login]])
            for login in missing(rc, list(types)) if types[login]]


def describe(types):
    return " and ".join(", ".join(f"{TYPES[t][0]} {t}" for t in types).rsplit(", ", 1))


# Pending credits are kept in the body of the open credit pull request, one marker
# per credit, so each new command rebuilds the branch from develop without conflicts.
MARKER = re.compile(rf"<!-- credit: (?P<login>{LOGIN}) (?P<types>[a-z,]+) #(?P<thread>\d+) -->")


def markers(body):
    """The pending credits in a pull request body, as (login, types, thread)."""
    return [(m["login"], m["types"].split(","), int(m["thread"])) for m in MARKER.finditer(body or "")]


def pending_body(credits):
    lines = ["Adds these people to the contributor list in the README:", ""]
    for login, types, thread in credits:
        lines.append(f"- @{login} for {describe(types)}, requested in #{thread}")
    lines += ["", "When this is merged, each person is thanked on the issue or pull request "
              "where they were credited. To add more people before merging, comment "
              "`/credit @user <type>` on any issue or pull request.", ""]
    lines += [f"<!-- credit: {login} {','.join(types)} #{thread} -->" for login, types, thread in credits]
    return "\n".join(lines) + "\n"


def apply(rc, body, comment, thread, user=github_user):
    """Apply the pending credits in `body` and the commands in `comment`.

    Returns (credits, new, errors): all pending credits as (login, types, thread),
    the credits this comment added as (login, new types), and the errors.
    """
    pending = markers(body)
    for login, types, _ in pending:
        add(rc, login, types, user(login) if not find(rc, login) else None)
    commands, errors = parse(comment)
    new = []
    for login, types in commands:
        try:
            info = None if find(rc, login) else user(login)
        except Exception:  # noqa: BLE001 - any failure means the login is unusable
            errors.append(f"`@{login}`: no such GitHub user")
            continue
        added = add(rc, login, types, info)
        login = find(rc, login)["login"]
        if added:
            new.append((login, added))
            pending.append((login, added, thread))
        else:
            errors.append(f"`@{login}` is already credited for {describe(types)}")
    return pending, new, errors


def thanks(body):
    """One thank-you comment per thread for a merged credit pull request."""
    by_thread = {}
    for login, types, thread in markers(body):
        by_thread.setdefault(thread, []).append((login, types))
    link = f"https://github.com/{REPO}#contributors"
    return {thread: "\n".join(
        f"@{login}, thank you! You are now in the [Human-GEM contributor list]({link}) "
        f"for {describe(types)}." for login, types in people)
        for thread, people in by_thread.items()}


def suggestion(credits):
    """The comment proposing credits after a merge, or "" if there are none."""
    if not credits:
        return ""
    lines = "\n".join(f"/credit @{login} {', '.join(types)}" for login, types in credits)
    return ("These people are not in the contributor list yet. To credit them, a maintainer "
            "replies with the lines below, changing the types if needed "
            f"({', '.join(f'`{t}`' for t in TYPES)}):\n\n```\n{lines}\n```")


def gh_json(*args):
    import subprocess
    out = subprocess.run(["gh", *args], capture_output=True, text=True, encoding="utf-8", check=True)
    return json.loads(out.stdout)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("add", help="credit a person and update the README")
    a.add_argument("login")
    a.add_argument("types", help="comma-separated, from: " + ", ".join(TYPES))
    sub.add_parser("generate", help="rebuild the README table from .all-contributorsrc")
    sub.add_parser("check", help="fail if the list has errors or the README is out of date")
    # Used by .github/workflows/add-contributor.yml. They read their input from the
    # environment, so that comment text never passes through the shell.
    w = sub.add_parser("apply", help="(workflow) apply COMMENT_BODY on top of PENDING_BODY")
    w.add_argument("--thread", type=int, required=True)
    w.add_argument("--out", required=True, help="JSON file for the result")
    s = sub.add_parser("suggest", help="(workflow) print credits to propose for a merged PR")
    s.add_argument("pr", type=int)
    sub.add_parser("thank", help="(workflow) print thank-you comments for PENDING_BODY as JSON")
    args = p.parse_args(argv)

    rc = load()
    if args.cmd == "add":
        types, unknown = parse_types(args.types)
        if unknown or not types:
            p.error(f"unknown type(s) {unknown}; use {', '.join(TYPES)}")
        new = add(rc, args.login, types)
        save(rc)
        generate(rc)
        print(f"{args.login}: {', '.join(new) if new else 'already credited for ' + ', '.join(types)}")
    elif args.cmd == "generate":
        generate(rc)
    elif args.cmd == "check":
        errors = problems(rc, README.read_bytes().decode("utf-8"))
        for e in errors:
            print(e, file=sys.stderr)
        return 1 if errors else 0
    elif args.cmd == "apply":
        pending, new, errors = apply(rc, os.environ.get("PENDING_BODY", ""),
                                     os.environ.get("COMMENT_BODY", ""), args.thread)
        if new:
            save(rc)
            generate(rc)
        Path(args.out).write_text(json.dumps({
            "new": [{"login": login, "types": types, "text": f"@{login} for {describe(types)}"}
                    for login, types in new],
            "errors": errors,
            "body": pending_body(pending),
            "title": "docs: credit " + ", ".join(sorted({login for login, _, _ in pending})),
        }, ensure_ascii=False), encoding="utf-8")
    elif args.cmd == "suggest":
        pr = gh_json("pr", "view", str(args.pr), "--repo", REPO,
                     "--json", "author,files,closingIssuesReferences")
        issues = []
        for ref in pr.get("closingIssuesReferences") or []:
            i = gh_json("issue", "view", str(ref["number"]), "--repo", REPO, "--json", "author,labels")
            issues.append((i["author"]["login"], [lab["name"] for lab in i["labels"]]))
        # People already waiting in the open credit pull request are not proposed again.
        for login, types, _ in markers(os.environ.get("PENDING_BODY", "")):
            add(rc, login, types, {"login": login, "avatar_url": "", "html_url": ""})
        print(suggestion(suggest(rc, pr["author"]["login"], [f["path"] for f in pr["files"]], issues)))
    elif args.cmd == "thank":
        print(json.dumps({str(k): v for k, v in thanks(os.environ.get("PENDING_BODY", "")).items()},
                         ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
