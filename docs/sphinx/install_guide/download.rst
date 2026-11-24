.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _download:


Download Siconos
================

Latest source release
---------------------

The latest source release can be downloaded in zip or tar.gz format from:

https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/releases
or
https://github.com/siconos/siconos/releases/latest


Development sources
-------------------

The Siconos project is hosted at https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos.git

and the development source code can be freely downloaded, for example with::

  git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos.git

As a user, you will probably only need to clone the repository (as shown above) once, and then update your local copy to
get the latest revision::

  cd /path/to/siconos
  git pull

As a developer, you need access to the private (devel) part of Siconos. Contact us if required.
You will also need to learn more about git, see for example https://git-scm.com/book/en/v1/Getting-Started-About-Version-Control.

Below is a short git refresher:

* bring your working copy up to date with the remote repository::

    git pull --rebase

* commit the new version of your file(s) to your local repository::

    git commit -a -m "some comments"

* check the status of your local repository::

    git status

* add a file to the index::

    git add filename

* remove a file from the index::

    git rm filename

* see the diff between your branch (here main) and another one (here the remote origin)::

    git diff origin main

* see the list of files that differ::

    git diff origin main --stat

* propagate your changes to the main repository::

    git push

Packages
--------

Guix
""""

https://guix.gnu.org/

Siconos is available as a Guix package. Please contact us for details about the required channels and how to install Siconos with the guix command.


Debian
""""""

.. warning::
   The Siconos version in the Debian packages is outdated.


* Siconos 4.2.0 https://packages.debian.org/sid/siconos  sid (unstable)
