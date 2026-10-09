######################
Documenting PyCBC code
######################

These page contains some details of how to use the pyCBC Sphinx documentation system. This is a very similar to that used by numpy and scipy to document that code.

All code in PyCBC must be appropriately documented. If you find any documentation to be inadequate, please contact the development team, or improve yourself where possible.

Sphinx
======

PyCBC uses the `sphinx` documentation system to produce automatic documentation from all modules and functions, with a layer on top to describe the code appropriately.

The `sphinx website <http://sphinx-doc.org>`_ has a `nice tutorial <http://sphinx-doc.org/rest.html#rst-primer>`_ on the reStructuredText format used for documentation.

Documenting the PyCBC package
=============================

The overview documentation pages (e.g. this page), live in the ``docs`` directory of the git repository. These pages are linked together using internal links.

For example: the file ``docs/index.rst`` file is used to build the home page for the documentation, and contains the following block::

    .. toctree::
        :maxdepth: 1

        install
        documentation

The ``install`` and ``documentation`` entries refer to the ``install.rst`` and ``documentation.rst`` files in the same directory.

In order to provide readable documentation for all users, each key module in ``pycbc`` should have its own directory in ``docs`` that builds upon the automatic documentation built from its classes and functions. For example, the ``pycbc.frame`` module has it's own documentation in ``docs/frame.rst`` that includes some custom text, and a copy of the automatic documentation from its classes and functions.

How much text goes directly into the module docstrings, and how much is abstracted into a separated ``rst`` file in the ``docs/`` directory is a matter of personal taste, and keeping the code readable, but should have little effect on the HTML documentation.

Documenting PyCBC modules
=========================

All PyCBC modules should be documented using ``sphinx`` syntax, in the same way that ``matplotlib``, ``numpy`` and ``scipy`` functions are listed.

The numpy github repository includes a `nice style guide <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_ for documenting modules and their members using sphinx.


Documenting PyCBC scripts
=========================

Documenting command-line scripts isn't ideal in any documentation language, including sphinx.

In ``PyCBC``, command-line scripts in the ``bin`` directory of the git repository should be accompanied by a reStructuredText file in ``docs/bin``.

However, at some point this directory got removed. So what is the current recommendation for documenting scripts??

This file contains an overview of what the code does, and some other information, with finally a dump of the command-line ``--help`` message via this directive::

    .. command-output:: pycbc_inspiral --help
