==========================
Running PyCBC under Docker
==========================

The easiest way to start using PyCBC is to install one of our `Docker containers <https://hub.docker.com/u/pycbc/>`_. First, install the `Docker Community Edition <https://www.docker.com/community-edition>`_ for your `Mac <https://store.docker.com/editions/community/docker-ce-desktop-mac?tab=description>`_ or `Windows <https://store.docker.com/editions/community/docker-ce-desktop-windows?tab=description>`_ desktop. Docker CE installations for `Linux platforms <https://www.docker.com/community-edition#/download>`_ are also available.


To start a Docker container with no graphics, type the commands::

    docker pull pycbc/pycbc-el7:latest
    docker run -it pycbc/pycbc-el7:latest

This example downloads current version of the code from the `GitHub master branch. <https://github.com/ligo-cbc/pycbc>`_ Replace the string ``latest`` with one of the `PyCBC release tags <https://github.com/ligo-cbc/pycbc/releases>`_ (e.g. ``v1.7.0``) to install a container containing a released version of PyCBC. The container includes all of the required software and dependencies to run PyCBC, including a compatible version of LALSuite installed into the root filesystem. The command above starts a login shell as the pycbc user. To override this and log in as root, run the command::

   docker run -it pycbc/pycbc-el7:latest /bin/bash -l

-------------------------------------
Using jupyter notebook within docker
-------------------------------------

One can start a jupyter notebook within docker and then port forward to your
computer's environment.::

    docker run -it -p 8888:8888 --name pycbc_test pycbc/pycbc-el7:latest /bin/su -l pycbc -c "jupyter notebook --no-browser --ip 0.0.0.0"

Once the image is running, you can connect from your computer's web browser to the address printed to the screen by jupyter. This is typically the local host adddress, e.g. ``127.0.0.1``

-------------------------------
Sharing user files and SSH keys
-------------------------------

It can be useful to share your SSH public/private key with the Docker container, for example to allow you to git push and pull from your repository on GitHub. To do this, add the argument ``-v ${HOME}/.ssh:/opt/pycbc/.ssh`` to the ``docker run`` commands.  You can also create e.g. a ``scratch`` directory and use the ``-v`` option to mount it in the container. This directory can be used to transfer files between the container and the host computer. See the `Docker volumes documentation <https://docs.docker.com/storage/volumes/>`_ for a detailed explaination of mounting directories inside a docker container.
