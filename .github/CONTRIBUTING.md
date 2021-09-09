## Contributing guidelines

Contributions to **Human-GEM** are very welcome and greatly appreciated! Credit will always be given to anyone who contribute.

You can contribute in **2** major ways: by creating issues, and by sending pull requests (PRs) with additions, deletions, corrections, etc. to the `yml` model and/or `tsv` annotation files. Please follow the following guidelines:

### Reporting issues

Report an issue at [here](https://github.com/SysBioChalmers/Human-GEM/issues), if you notice any of the following:

* Incorrect annotation for any model components.
* Missing feature or field you would like the model to have.
* Bug/weird simulation results.
* Lacking documentation.
* Any type of feedback.

If you are unsure about the issue, consider asking first in our [Gitter chat room](https://gitter.im/SysBioChalmers/Human-GEM).

When creating the issue, please make sure:

* You checked that a similar issue does not exist already
* You tested your code (if any) with all requirements for running the model.
* You did your analysis in the `main` branch of the repository.
* You provide any necessary files/links needed for understanding the issue.

Feel free to also comment on any of the [open issues](https://github.com/SysBioChalmers/Human-GEM/issues).


### Contributing to the model/annotation

If you want to contribute to the model with some additions or improvements, please always start by raising an issue for asking/describing what you want to achieve. This way, we reduce the risk of duplicated efforts and you may also get suggestions on how to best proceed, e.g. there may be half-finished work in some branch that you could work with. Also, feel free to browse our [open issues](https://github.com/SysBioChalmers/Human-GEM/issues) and our [ongoing projects](https://github.com/SysBioChalmers/Human-GEM/projects): Anything tagged with "help wanted" is open to whoever wants to implement it.


Here's how to set up Human-GEM to contribute small features or changes. This should work for most curation work.

1. Fork the [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) repository to your local, by clicking on the upper right corner. Then, switch to the `develop` branch before starting to work in the forked repository.

2. Modify model file `Human-GEM.yml` and/or annotation files `reactions.tsv`, `metabolites.tsv`, and `genes.tsv`, and commit changes as suggested below in "Semantic commits". 

3. Submit a pull request from the `develop` branch of the forked repo on GitHub website to the `develop` branch of the original repo. We recommend ticking the box "Allow edits from maintainers" if you wish for us to be able to contribute directly to your branch (speeding-up the reviewing process).


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


## Acknowledgments

These contribution guidelines were adapted from the guidelines of [yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM/blob/main/.github/CONTRIBUTING.md).

