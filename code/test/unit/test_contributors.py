"""Tests for code/contributors.py, which keeps the contributor list and README table."""

import contributors as c


def person(login, name=None, types=("bug",)):
    return {"login": login, "name": name or login, "avatar_url": f"https://a/{login}",
            "profile": f"https://github.com/{login}", "contributions": list(types)}


def rc(*people):
    return {"contributorsPerLine": 2, "badgeTemplate": "[badge-<%= contributors.length %>]",
            "contributors": list(people)}


def user(login):
    return {"login": login, "name": None, "avatar_url": f"https://a/{login}",
            "html_url": f"https://github.com/{login}", "blog": ""}


def test_parse_credit_command_and_legacy_phrase():
    body = ("Thanks!\n/credit @ann bug, research\n/credit bob for code and docs\n"
            "@all-contributors please add @cy for ideas.")
    credits, errors = c.parse(body)
    assert credits == [("ann", ["bug", "research"]), ("bob", ["code", "doc"]), ("cy", ["ideas"])]
    assert errors == []


def test_parse_reports_unknown_or_missing_types():
    credits, errors = c.parse("/credit @ann content\n/credit @bob")
    assert credits == []
    assert "`content`" in errors[0]
    assert len(errors) == 2


def test_parse_ignores_quoted_and_inline_mentions():
    assert c.parse("> /credit @ann bug\nplease /credit @bob bug") == ([], [])


def test_add_merges_types_in_fixed_order():
    r = rc(person("ann", types=["code"]))
    assert c.add(r, "ANN", ["bug", "code"]) == ["bug"]
    assert r["contributors"][0]["contributions"] == ["bug", "code"]


def test_add_new_person_uses_github_profile():
    r = rc()
    c.add(r, "bob", ["review"], user("bob"))
    assert r["contributors"] == [{"login": "bob", "name": "bob", "avatar_url": "https://a/bob",
                                  "profile": "https://github.com/bob", "contributions": ["review"]}]


def test_render_sorts_by_name_and_updates_badge():
    readme = f"[badge-0]x\n{c.LIST_START}\nold\n{c.LIST_END}\n"
    readme = readme.replace("[badge-0]", "<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify "
                            "this section -->\n[badge-0]\n<!-- ALL-CONTRIBUTORS-BADGE:END -->")
    out = c.render(rc(person("zed", "Zed"), person("amy", "Amy"), person("bo", "Bo")), readme)
    assert "[badge-3]" in out and "old" not in out
    assert out.index('alt="Amy"') < out.index('alt="Bo"') < out.index('alt="Zed"')
    assert out.count("<tr>") == 2
    assert "issues?q=author%3Aamy" in out


def test_problems_finds_duplicates_and_unknown_types():
    r = rc(person("ann"), person("Ann", types=["content"]))
    found = c.problems(r, "")
    assert any("more than once" in p for p in found)
    assert any("unknown types" in p for p in found)


def test_repository_list_is_valid_and_readme_in_sync():
    assert c.problems(c.load(), c.README.read_bytes().decode("utf-8")) == []


def test_suggest_only_people_not_listed():
    r = rc(person("maintainer"))
    got = c.suggest(r, "newdev", ["model/reactions.tsv", "code/x.py", "README.md"],
                    [("reporter", ["bug"]), ("asker", ["feature"]), ("maintainer", ["bug"]),
                     ("dependabot[bot]", [])])
    assert got == [("reporter", ["bug"]), ("asker", ["ideas"]),
                   ("newdev", ["research", "code", "doc"])]
    assert c.suggest(r, "maintainer", ["code/x.py"], []) == []


def test_suggestion_text_has_copyable_lines():
    text = c.suggestion([("ann", ["bug", "research"])])
    assert "```\n/credit @ann bug, research\n```" in text
    assert c.suggestion([]) == ""


def test_apply_rebuilds_pending_and_adds_new():
    r = rc(person("old"))
    body = c.pending_body([("ann", ["bug"], 5)])
    pending, new, errors = c.apply(r, body, "/credit @bob code\n/credit @old bug", 9, user=user)
    assert pending == [("ann", ["bug"], 5), ("bob", ["code"], 9)]
    assert new == [("bob", ["code"])]
    assert "already credited" in errors[0]
    assert {p["login"] for p in r["contributors"]} == {"old", "ann", "bob"}
    assert c.markers(c.pending_body(pending)) == pending


def test_apply_reports_unknown_user():
    def missing_user(login):
        raise OSError("404")
    _, new, errors = c.apply(rc(), "", "/credit @ghost bug", 1, user=missing_user)
    assert new == [] and "no such GitHub user" in errors[0]


def test_thanks_one_comment_per_thread():
    body = c.pending_body([("ann", ["bug"], 5), ("bob", ["code", "doc"], 5), ("cy", ["ideas"], 7)])
    t = c.thanks(body)
    assert set(t) == {5, 7}
    assert "@ann" in t[5] and "@bob" in t[5] and "💻 code and 📖 doc" in t[5]
