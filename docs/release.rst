##########################
Creating Releases of PyCBC
##########################

.. note::

    Only authorized maintainers of PyCBC should use these instructions.


==============================
Creating the release on GitHub
==============================

To create a new PyCBC release:

#. Make sure that the setup.py file contains the correct version number, which should be in the format ``x.y.x`` (where x, y, and z are the major, minor, and patch levels of the release) in PyCBC's setup.py file.
#. Set ``release = True`` in the PyCBC's setup.py file.
#. Commit the changed setup.py file and push to commit to the repository.
#. Go to the `PyCBC release page <https://github.com/ligo-cbc/pycbc/releases>`_ and click on ``Draft a new release``. 
#. Enter a tag version in the format ``vx.y.z``. Note the ``v`` in front of the major, minor, and patch numbers. 
#. Enter a descriptive release title and write a desciption of the release in the text box provided.
#. Click on ``Publish release`` to create the release.
#. Update the setup.py file with an incremented major or minor version number and make sure that the string ``dev`` appears in that version. For example, if you just released ``1.2.1`` then change the string to ``1.3.dev0`` or ``2.0.dev0`` as appropriate. This is need to ensure that if someone is building from source, it always takes precedence over an older release version.
#. Set ``release = False`` in PyCBC's setup.py file.
#. Commit these changes and push them to the repository.

=====================================
Uploading to the Python Package Index
=====================================

.. note::

    Keep in mind that once you register a version number with the package index and upload a package, you can never change it. You can only deactivate the package from the index, and increment the version number.

If this is the first time that you have pushed to PyPI, you will need to create a configuration file. Create a file in your home directory called ``.pypirc`` that contains the lines

.. highlight:: 

    [distutils] # this tells distutils what package indexes you can push to
    index-servers =
      pypi
      pypitest

    [pypi]
    repository: https://pypi.python.org/pypi
    username: your_username
    password: your_password

    [pypitest]
    repository: https://testpypi.python.org/pypi
    username: your_username
    password: your_password

Replace ``your_username`` and ``your_password`` with your PyPI and PyPI testing usernames and passwords.

Once you have tagged the release of PyCBC a tarball should be available from the `GitHib PyCBC releases page <https://github.com/ligo-cbc/pycbc/releases>`_. Download this tarball, untar it, and change into the source directory of the tarball. 

First check the release with the PyPI test repository. Register the package with 

.. code:: bash

    python setup.py register -r pypitest

If you get no errors, you can then upload the package with

.. code:: bash

    python setup.py sdist upload -r pypitest

You should now see the package uploaded in the `PyPI Test Repository <https://testpypi.python.org/pypi>`_. If this is successful, you can puhlish it to the main repository with the commands

.. code:: bash

    python setup.py register -r pypi
    python setup.py sdist upload -r pypi

The package should then be available in PyPI for download.
