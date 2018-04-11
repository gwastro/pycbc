# Contributing to PyCBC

This page outlines the recommended procedure for contributing changes to the PyCBC repository. Please familiarise yourself with [GitHub](https://github.com) ensure your account is configured according to these instructions.

## Reporting Issues

When reporting issues, please include as much detail as possible to reproduce the error, including information about your operating system and the version of each (relevant) component of PyCBC.
If possible, please include a brief, self-contained code example that demonstrates the problem.

## Contributing code

All contributions to PyCBC code must be made using the [GitHub Flow](https://guides.github.com/introduction/flow/) model, which must then be reviewed by one of the project maintainers.

If you wish to contribute new code, or changes to existing code, please follow the following development workflow.

### Make a fork (copy) of PyCBC

**You only need to do this once**

1. Go to the [PyCBC repository home page](https://github.com/gwastro/pycbc)
2. Click on the *Fork* button (top-right-hand corner)
3. Select the namespace that you want to create the fork in, this will usually be your personal namespace

### Clone your fork

```bash
git clone https://github.com/<username>/pycbc.git
```

### Updating your fork

If you already have a fork of PyCBC, and are starting work on a new project you can link your clone to the main (`gwastro`) repository and pull in changes that have been merged since the time you created your fork, or last updated:

1. Link your fork to the main repository:

    ```bash
    cd pycbc
    git remote add gwastro https://github.com/gwastro/pycbc.git
    ```

2. Fetch new changes from the `gwastro` repo

    ```bash
    git fetch gwastro
    ```

### Creating a new feature branch

All changes should be developed on a feature branch, in order to keep them separate from other work, simplifying review and merge once the work is done.

To create a new feature branch:

```bash
git fetch gwastro
git checkout -b my-new-feature gwastro/master
```

### Hack away

1. Develop the changes you would like to introduce, using `git commit` to finalise a specific change.
   Ideally commit small units of change often, rather than creating one large commit at the end, this will simplify review and make modifying any changes easier.

    Commit messages should be clear, identifying which code was changed, and why.
   Common practice is to use a short summary line (<50 characters), followed by a blank line, then more information in longer lines.

2. Push your changes to the remote copy of your fork on GitHub

    ```bash
    git push origin my-new-feature
    ```
   **Note:** For the first `push` of any new feature branch, you will likely have to use the `-u/--set-upstream` option to `push` to create a link between your new branch and the `origin` remote:

    ```bash
    git push --set-upstream origin my-new-feature
    ```

### Open a Pull Request

When you feel that your work is finished, you should create a Pull Request to propose that your changes be merged into the main (`gwastro`) repository.

After you have pushed your new feature branch to `origin`, you should find a new button on the [PyCBC repository home page](https://github.com/gwastro/pycbc/) inviting you to create a Pull Request out of your newly pushed branch.
You should click the button, and proceed to fill in the title and description boxes on the PR page.

Once the request has been opened, one of the maintainers will assign someone to review the change.

## More Information

More information regarding the usage of GitHub can be found in the [GitHub Guides](https://guides.github.com/).
