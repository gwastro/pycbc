##########################
Creating Releases of PyCBC
##########################

..note::

    Only authorized maintainers of PyCBC should use these instructions.


==============================
Creating the release on GitHub
==============================

To create a new PyCBC release:

    1. Make sure that the setup.py file contains the correct version number, which should be in the format ``x.y.x`` (where x, y, and z are the major, minor, and patch levels of the release) in PyCBC's setup.py file.
    1. Set ``release = True`` in the PyCBC's setup.py file.
    1. Commit the changed setup.py file and push to commit to the repository.
    1. Go to the `PyCBC release page <https://github.com/ligo-cbc/pycbc/releases>`_ and click on ``Draft a new release``. 
    1. Enter a tag version in the format ``vx.y.z``. Note the ``v`` in front of the major, minor, and patch numbers. 
    1. Enter a descriptive release title and write a desciption of the release in the text box provided.
    1. Click on ``Publish release`` to create the release.
    1. Update the setup.py file with an incremented major or minor version number and make sure that the string ``dev`` appears in that version. For example, if you just released ``1.2.1`` then change the string to ``1.3.dev0`` or ``2.0.dev0`` as appropriate. This is need to ensure that if someone is building from source, it always takes precedence over an older release version.
    1. Set ``release = False`` in PyCBC's setup.py file.
    1. Commit these changes and push them to the repository.

=====================================
Uploading to the Python Package Index
=====================================

Follow the `instructions for submitting a package to PyPI <http://peterdowns.com/posts/first-time-with-pypi.html>`_ to submit the package to the Python Package Index. 

Keep in mind that once you register a version number with the package index and upload a package, you can never change it. You can only deactivate the package from the index, and increment the version number.

