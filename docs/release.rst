##########################
Creating Releases of PyCBC
##########################

.. note::

    Only authorized maintainers of PyCBC should use these instructions.


==============================
Creating the release on GitHub
==============================

To create a new PyCBC release:

#. Go to the `PyCBC release page <https://github.com/gwastro/pycbc/releases>`_ and click on ``Draft a new release``.
#. Enter a tag version in the format ``vx.y.z``. Note the ``v`` in front of the major, minor, and patch numbers.
#. Enter a descriptive release title and write a description of the release in the text box provided.
#. Click on ``Publish release`` to create the release.

.. note::

    Releases of PyCBC should be made from the master. Do not make a branch
    unless you are back-porting a bug fix from a new release series to an
    old production release series.

Creating the release will trigger a CI build that updates CVMFS, Docker, and PyPI with the release, and increments the patch number.
Please ensure that you check the outputs of these builds and urgently report any errors to other maintainers, you will be emailed if the builds fails.

------------------------------------------------
Backporting Bug Fixes to Previous Release Series
------------------------------------------------

Branches should only be created when bug fixes from master need to be back
ported to an old release series (e.g. adding a bug fix from the 1.4 series to
the 1.3 series. )

To create a branch for the bug fix, make a new branch from the last release
tag and cherry pick the changes to that branch. For example, to create a
``1.3.8`` branch to add bug fixes to the ``1.3.7`` release, you should

.. code:: bash

  git checkout -b b1.3.8 v1.3.7
  git cherry-pick [hash of first change]
  git cherry-pick [hash of second change]
  [continue using git cherry-pick for all new changes]
  git commit
  git push

Then go to the `PyCBC release page <https://github.com/gwastro/pycbc/releases>`_ and click on ``Draft a new release`` following the instructions above.

