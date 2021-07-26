## Contributor guidelines

Contributions to Human-GEM are very welcome and greatly appreciated! Credit will always be given to anyone who contribute.

You can contribute in 2 main ways: by creating issues, and by sending pull requests (PRs) with additions, deletions, corrections, etc. to the model. Please follow the following guidelines:

### Reporting issues in the model

Report an issue at https://github.com/SysBioChalmers/Human-GEM/issues if you note any of the following:

* Incorrect annotation for any model component.
* Missing feature or field you would like the model to have.
* Bug/weird simulation results.
* Lacking documentation.
* Any type of feedback.

If you are unsure about the issue, consider asking first in our [Gitter chat room](https://gitter.im/SysBioChalmers/Human-GEM).

When creating the issue, please make sure:

* You tested your code (if any) with all requirements for running the model.
* You did your analysis in the `main` branch of the repository.
* You provide any necessary files/links needed for understanding the issue.
* You checked that a similar issue does not exist already

Feel free to also comment on any of the [open issues](https://github.com/SysBioChalmers/Human-GEM/issues).


### Contributing to the model

If you want to contribute to the model with some additions or improvements, please consider starting by raising an issue and assign it to yourself to describe what you want to achieve. This way, we reduce the risk of duplicated efforts and you may also get suggestions on how to best proceed, e.g. there may be half-finished work in some branch that you could start with. Also, feel free to browse our [open issues](https://github.com/SysBioChalmers/Human-GEM/issues) and our [ongoing projects](https://github.com/SysBioChalmers/Human-GEM/projects): Anything tagged with "help wanted" is open to whoever wants to implement it.

Here's how to set up Human-GEM for local development to contribute smaller features or changes that you can implement yourself:

1. Make sure that you have all [requirements](https://github.com/SysBioChalmers/Human-GEM#required-software) for contributing to Human-GEM. Note that COBRA and RAVEN should be cloned to your computer and not just downloaded.

2. Fork the Human-GEM repository on GitHub (go to https://github.com/SysBioChalmers/Human-GEM & click on the upper right corner).

3. Clone your fork locally:
    ```
    git clone https://github.com/<your Github name>/Human-GEM.git
    ```
	
4. Check out the branch that you want to contribute to. Most likely that will be `devel`:
    ```
    git checkout devel
    ```
	
5. Create a branch for local development based on the previously checked out branch ([see below](#branching-model) for details on the branching model and how to name your branch):
    ```
    git checkout -b name-of-your-branch
    ```
	
6. Now you can make your changes locally!
    * If your changes are minor (e.g. a single chemical formula you wish to correct), you can do it directly from the command line.
    * If your changes are not so small and require several steps, create a script that loads the model, reads data (if applicable), changes the model accordingly, and saves the model back.
    * Each script should start with a commented section describing the script, explaining the parameters, and indicating your name and the date it was written. Existing functions can clarify what style should be used.
    * Store scripts in the appropriate folder in `/code` and data (as `.tsv` files) in the appropriate folder in `/data`. If you think no folder is adequate for your script/data, feel free to create your own folder. Note that binary data such as `.mat` structures or `.xls` tables cannot be stored in the repo (as they cannot be version-controlled, and they increment too much the size of the repo).
    * When you are done making changes, review locally your changes with `git diff` or any git client, to make sure you are modifying the model as you intended.

7. Commit your changes and push your branch to GitHub.
    ```
    git add .
    git commit -m "Title of your commit"
    git push origin name-of-your-branch
    ```
    [See below](#semantic-commits) for recommendations on how to name your commits. In case of larger updates, you can of course make several commits on a single contribution. However, if you need to make many commits, consider if your contribution could be instead split into separate branches (making it easier for reviewing later).

8. Submit a pull request through the GitHub website (https://help.github.com/articles/creating-a-pull-request-from-a-fork/) to the `devel` branch of the original SysBioChalmers repo (not to your fork). We recommend ticking the box "Allow edits from maintainers" if you wish for us to be able to contribute directly to your branch (speeding-up the reviewing process).

Finally, and for larger features that you want to work on collaboratively, you may consider to first request to join our development team to get write access to the repository so that you can create a branch directly in the main repository (or simply ask the administrator to create a branch for you). Once you have a new branch, you can push your changes directly to the main repository and when finished, submit a pull request from that branch to `devel`. [See below](#development-team-guidelines) for more details.

Thank you very much for contributing to Human-GEM!

#### Branching model

* `devel`: Is the branch to which all pull requests should be made.

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

These contribution guidelines were adapted from the guidelines of [SysBioChalmers/yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM/blob/main/.github/CONTRIBUTING.md).
