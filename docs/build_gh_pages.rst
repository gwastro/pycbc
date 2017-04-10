########################################
Building Documentation for Git Hub Pages
########################################

===============================
Creating a Git Hub project page
===============================

Git Hub pages are built from a branch of the repository called ``gh-pages``.
If you have not already created a Git Hub project pages for PyCBC in your
repository, make a ``gh-pages`` branch in your repository as follows:

.. code-block:: bash

    git checkout --orphan gh-pages
    git rm -rf .
    git clean -dxf
    touch .nojekyll
    mkdir latest
    git add .nojekyll latest
    git commit -a -m "set up gh-pages branch"
    git push origin gh-pages
    git branch --set-upstream gh-pages origin/gh-pages

These commands create the branch and then remove all of the files from this branch, so that it will just contain the documentation pages.

.. note::

    The main `ligo-cbc/pycbc <https://github.com/ligo-cbc/pycbc>`_ repository already has a `gh-pages` branch, so do not do this in the main repository.

======================================
Building and pushing the documentation
======================================

The documentation should built from the source code on a regular branch and installed into the ``gh-pages`` branch. Since git cannot have two branches checked out simultaneously, you need to check out another copy of the repository with the ``gh-pages`` branch inside your PyCBC source repository. Do this will the following commands (assuming your PyCBC git repository is in a directory named ``pycbc``).

.. code-block:: bash

    cd /path/to/your/repo/pycbc
    git clone git@github.com:github-username/pycbc.git _gh-pages

Now flush the contents of this directory. We do this, as the documentation is
not really under version control in the ``gh-pages`` branch. We just use this
branch to publish to GitHub pages. Run the commands

.. code-block:: bash

    cd _gh-pages
    git rm -rf *
    git commit -a -m "flush documentation"
    cd ..

The last ``cd`` command should put you back at the top level of your PyCBC
source directory. To build the documentation into the ``_gh-pages`` directory,
run the command

.. code-block:: bash

    python setup.py build_gh_pages
    
This will build the documentation in the second repository that you created called ``_gh-pages`` under the directory ``latest/``. To push these changes up to Git Hub

.. code-block:: bash

    cd _gh-pages
    git add --all
    git commit -a -m "documentation update"
    git push origin gh-pages

The documentation will then be available under your Git Hub pages at ``http://username.github.io/pycbc/latest/html/`` where you should replace ``username/`` with your Git Hub account name.

In this example, we checkout the master branch to build the documentation, but you can change the last command above to checkout any other branch that your are developing.

.. note::

    Be careful with the ``git rm -rf *`` command as if you run it in the wrong
    directory you can delete the contents of your git repository. If you do
    this by accident, you can use ``git reset`` to undo the commit.

