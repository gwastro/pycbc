######################
Documenting PyCBC code
######################

All code in PyCBC should be appropriately documented.

Sphinx
======

PyCBC uses the `sphinx` documentation system to produce automatic documentation from all modules and functions, with a layer on top to describe the code appropriately.

Documenting PyCBC modules
=========================

All PyCBC modules should be documented using `sphinx` syntax, in the same way that `matplotlib`, `numpy` and `scipy` functions are listed.

Documenting PyCBC scripts
=========================

Documenting command-line scripts isn't ideal in any documentation language, including `sphinx`. Command-line scripts in the `bin` directory of the git repository should be accompanied by a reStructuredText file in `docs/bin`.

For example, the `pycbc_inspiral` script is accompanied by `pycbc_inspiral.rst` which documents what the code does, and how to run it.

Example: :doc:`bin/pycbc_inspiral`
