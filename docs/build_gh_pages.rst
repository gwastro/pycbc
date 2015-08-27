################
Building Documentation for Git Hub Pages
################

=============
Creating a Git Hub project page
=============

Git Hub pages are built from a branch of the repository called ``gh-pages``.
If you have not already created a Git Hub project pages for PyCBC in your
repository, make a ``gh-pages`` branch in your repository as follows:

.. code-block:: bash

    git checkout --orphan gh-pages
    git rm -rf .
    git clean -dxf
    touch .nojeckyl
    mkdir latest
    git add .nojeckyl latest
    git commit -a -m "set up gh-pages branch"
    git push origin gh-pages

These commands create the branch and then remove all of the files from this branch, so that it will just contain the documentation pages.

.. note::

    The main `ligo-cbc/pycbc <https://github.com/ligo-cbc/pycbc>`_ repository already has a `gh-pages` branch, so do not do this in the main repository.

=============
Building and pushing the documentation
=============

The documentation should built from the source code on a regular branch and installed into the ``gh-pages`` branch. Since git cannot have two branches checked out simultaneously, you need to make a copy of your repository with just the gh-pages branch checked out. Do this will the following commands (assuming your PyCBC git repository is in a directory named ``pycbc``).

.. code-block:: bash

    cd /path/to/your/repo/pycbc
    git checkout gh-pages
    cd ..
    cp -a pycbc pycbc-gh-pages
    cd pycbc
    git checkout master

In this example, we checkout the master branch to build the documentation, but you can change the last command above to checkout any other branch that your are developing.  To generate the documentation, from the top level of the PyCBC source tree run

.. code-block:: bash

    python setup.py build_gh_pages
    
This will build the documentation in the second repository that you created called ``pycbc-gh-pages`` under the directory ``latest/``. To push these changes up to Git Hub

.. code-block:: bash

    cd ../pycbc-gh-pages
    git add --all
    git commit -a -m "updated documentation"
    git push

The documentation will then be available under your Git Hub pages at ``http://username.github.io/pycbc/latest/`` where you should replace ``username/`` with your Git Hub account name.
