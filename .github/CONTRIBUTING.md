## Contributing guidelines

Contributions to **Human-GEM** are very welcome and greatly appreciated! Credit is given to everyone who contributes, tracked with [all-contributors](https://allcontributors.org); see [Crediting contributors](#crediting-contributors) below.

You can contribute in **2** major ways: by creating issues, and by sending pull requests (PRs) with additions, deletions, corrections, etc. to the `yml` model and/or `tsv` annotation files. Please follow the following guidelines:

### Reporting issues

Report an issue at [here](https://github.com/SysBioChalmers/Human-GEM/issues), if you notice any of the following:

* Incorrect annotation for any model components.
* Missing feature or field you would like the model to have.
* Bug/weird simulation results.
* Lacking documentation.
* Any type of feedback.

When creating the issue, please make sure:

* You checked that a similar issue does not exist already
* You tested your code (if any) with all requirements for running the model.
* You did your analysis in the `main` branch of the repository.
* You provide any necessary files/links needed for understanding the issue.

Feel free to also comment on any of the [open issues](https://github.com/SysBioChalmers/Human-GEM/issues).

To propose a specific correction or addition to the model, use the **Curation proposal** issue template. It asks for the current state, the evidence, the proposed change, counter-arguments and the expected effect.

### Labels and issue lifecycle

| Label | Meaning |
| --- | --- |
| `bug` | The model or code is wrong. |
| `feature` | Something new: a reaction, pathway, dataset or function. |
| `discussion` | A decision is needed before anyone implements a change. |
| `question` | A question about the model or its use. |
| `wip` | A pull request for this issue is open. |
| `fixed in develop` | The change is merged into `develop` and will be in the next release. |

An issue moves through these steps:

1. A maintainer sets `bug`, `feature`, `question` or `discussion`.
2. When a pull request with `Closes #nnn` in its description is opened against `develop`, issue `#nnn` gets `wip`.
3. When that pull request is merged, `wip` is replaced by `fixed in develop`. The issue stays open, because the change is not yet released.
4. At the next release, `fixed in develop` is removed and the issue is closed.

Steps 2 and 3 are done by the `Issue labels` workflow (`.github/workflows/issue-labels.yml`), and step 4 by the `Release` workflow (`.github/workflows/release.yml`), which comments on each closed issue with the release that contains its fix. An issue that will not be implemented is closed as "not planned", and a duplicate is closed as "duplicate" with a link to the other issue.


### Contributing to the model/annotation

If you want to contribute to the model with some additions or improvements, please always start by raising an issue for asking/describing what you want to achieve (for model changes, use the Curation proposal template). This way, we reduce the risk of duplicated efforts and you may also get suggestions on how to best proceed, e.g. there may be half-finished work in some branch that you could work with. Also, feel free to browse our [open issues](https://github.com/SysBioChalmers/Human-GEM/issues) and our [ongoing projects](https://github.com/SysBioChalmers/Human-GEM/projects): Anything tagged with "help wanted" is open to whoever wants to implement it.


Here's how to set up Human-GEM to contribute small features or changes. This should work for most curation work.

1. Fork the [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) repository to your local, by clicking on the upper right corner. Then, switch to the `develop` branch before starting to work in the forked repository.

2. Modify model file `Human-GEM.yml` and/or annotation files `reactions.tsv`, `metabolites.tsv`, and `genes.tsv`, and commit changes as suggested below in "Semantic commits". 

3. Submit a pull request from the `develop` branch of the forked repository to the `develop` branch of the original `SysBioChalmers/Human-GEM` repository. We recommend ticking the box "Allow edits from maintainers" if you wish for us to be able to contribute directly to your branch (speeding-up the reviewing process).


Finally, and for larger features that you want to work on collaboratively, you may consider to first request to join our development team to get write access to the repository so that you can create a branch directly in the main repository (or simply ask the administrator to create a branch for you). Once you have a new branch, you can push your changes directly to the main repository and when finished, submit a pull request from that branch to `develop`. [See below](#development-team-guidelines) for more details.

Thank you very much for contributing to Human-GEM!

#### Branching model

* `develop`: Is the branch to which all pull requests should be made.

* `main`: Is only modified by the administrator and is the branch with the tested & reviewed model that is released or ready for the next release.

* `{chore, doc, feat, fix, refactor, style}/descriptive-name`: Any other branch created in the model. If you work on a fix, start the branch name with `fix/`, if you work on a feature, start the branch name with `feat/`. Examples: `fix/format_reactions` or `feat/new_algorithms`. [See below](#semantic-commits) for more details on the possible actions you can use.
	
#### Semantic commits

Please use concise descriptive commit messages. Ideally, use semantic commit messages to make it easier to show what you are aiming to do:

`action: brief description`

`action` refers to what exactly are you doing in the commit, following a [standard definition](http://karma-runner.github.io/2.0/dev/git-commit-msg.html) in software development: 
* `chore`: updating toolbox, data files, etc.
* `doc`: updating documentation or explanatory comments in functions.
* `feat`: new feature added, e.g. new reaction / metabolite / function / etc.
* `fix`: something that was incorrect in the model and now has been corrected.
* `refactor`: see [code refactoring](https://en.wikipedia.org/wiki/Code_refactoring).
* `style`: minor format changes of model, functions or data (spaces, semi-colons, etc., no code change).

Some examples:

|commit|commit message|
|:---:|:---:|
|Add new rxns|`feat: add methanol pathway`|
|Remove a metabolite|`fix: remove duplicated citrate`|
|Add metabolite formula|`feat: add carbohydrate formulas`|
|Format name of compartment|`style: convert to all lowercase`|
|Split a rxn in 2|`refactor: split isomerase in 2 steps`|
|Update documentation of function|`doc: addDBnewRxn.m update documentation`|

More examples [here](https://github.com/SysBioChalmers/Human-GEM/commits/main). A more detailed explanation or comments is encouraged to be left in the commit description.


### Crediting contributors

Everyone who contributes is credited in the [Contributors](../README.md#contributors) section of the `README`, in the style of [all-contributors](https://allcontributors.org). That includes reporting an error in the model, not only writing code: most of the corrections in Human-GEM started as an issue.

**If you contributed**, you do not have to do anything. When your issue is fixed or your pull request is merged, a maintainer credits you, and you get a thank-you on that issue or pull request once you are in the list. If you think you were missed, say so in any issue.

**Maintainers** credit someone by commenting on any issue or pull request, one line per person:

```
/credit @username bug, research
```

Human-GEM uses six types:

| Type | | Used for |
|---|:--:|---|
| `bug` | 🐛 | reporting an incorrect reaction, metabolite, GPR, or annotation |
| `research` | 🔬 | evidence, values, and datasets behind a curation |
| `code` | 💻 | code, tests, and tooling |
| `review` | 👀 | reviewing pull requests |
| `doc` | 📖 | documentation and tutorials |
| `ideas` | 🤔 | ideas, planning, and coordination |

What happens then (`.github/workflows/add-contributor.yml`):

1. **After a merge into `develop`**, the workflow comments on the pull request with a ready-made `/credit` line for each person involved who is not in the list yet: the authors of the issues it closes (`bug`, or `ideas` for a `feature` issue) and the author of the pull request (`research` for `model/` or `data/`, `code` for code, `doc` for Markdown). Nothing is credited until a maintainer replies with those lines, edited where needed. People already in the list are not suggested again; add a new type for them with `/credit` yourself.
2. **On a `/credit` comment**, the person is added to a single open pull request, `contributors/pending`, against `develop`. The answer is added under the command: a link to that pull request, or what was wrong (an unknown type, an unknown user, already credited). Later `/credit` comments are added to the same pull request until it is merged, so there is one pull request to review rather than one per person.
3. **When that pull request is merged**, each person is thanked on the issue or pull request where they were credited, with a link to the list.

The list is kept in `.all-contributorsrc`, and the table in the `README` is generated from it by `code/contributors.py`. The `README` groups people in tiers by the number of contribution types they are credited for (`tiers` in `.all-contributorsrc`). To change the list by hand (for example a renamed GitHub account), edit `.all-contributorsrc` and run `python code/contributors.py generate`; `python code/contributors.py check` fails if a type is unknown, a person is listed twice, or the `README` is out of date. `python code/contributors.py add <username> <types>` does what `/credit` does, without the pull request.

The older `@all-contributors please add @username for bug` is still understood. The workflow acts on comments as it is on `main`, so a change to it takes effect after the next release.

### Releasing a new version

Releases are cut by the `Release` workflow (`.github/workflows/release.yml`), which runs the version bump, opens the release pull request, and publishes the release. The release notes are written by hand, in the description of the release pull request, and become the text of the GitHub release and of its Zenodo record; they are not kept as a file in the repository.

1. Start the workflow from the Actions tab, giving the new version (for example `2.1.0`). It must be a major, minor or patch increment of `version.txt` on `main`.
2. The workflow cuts `release/<version>` from `develop`, stamps the version, regenerates `model/Human-GEM.{yml,mat,xml,xlsx,txt}` and `model/dependencies.txt` with `code/io/increaseHumanGEMVersion.py`, fills the `{{nRXN}}` / `{{nMET}}` / `{{nGENE}}` placeholders in `README.md`, and opens a pull request into `main`. The Model QC workflow comments on that pull request as it does on any other, so the release is reviewed against the same checks.
3. Write the release notes in the description of that pull request, replacing its placeholder line. A check keeps the pull request from being merged while the notes are missing. See the wiki page [How to make a new release](https://github.com/SysBioChalmers/Human-GEM/wiki/How-to-make-a-new-release) for their structure, and earlier [releases](https://github.com/SysBioChalmers/Human-GEM/releases) for examples.

Merging that pull request tags `v<version>`, publishes the release titled `Human <version>` with the notes as its text, and opens a sync-back pull request into `develop`. The sync-back keeps the generated model files and `version.txt` off `develop`, and restores the `README.md` placeholders, so `develop` stays a template that is ahead of the last release rather than at a version.

## Acknowledgments

These contribution guidelines were adapted from the guidelines of [yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM/blob/main/.github/CONTRIBUTING.md).

