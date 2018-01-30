.. _download:

Download Siconos
================


Latest source release
---------------------

Latest source release is downloadable in zip or tar.gz format from here:

https://github.com/siconos/siconos/releases/latest


Development sources
-------------------

Siconos project is hosted on github : https://github.com/siconos/siconos

and the development source code can be freely downloaded. Try for example:

  git clone git@github.com:siconos/siconos.git

As user, you will probably only need to clone the repository (as shown above) once and then just update your local copy to
include the last revision::

  cd path-to-siconos
  git pull

As developer, you will need to learn more about git. Check for example https://git-scm.com/book/en/v1/Getting-Started-About-Version-Control.

Below, you can find a short git refresher:

* bring your working copy "up-to-date" with the github repository::

    git pull --rebase

* commit the new version of your file(s) to your local repository::

    git commit -a -m "some comments"

* check the status of your local repository::

    git status

* add a file to the index::

    git add filename

* remove a file from the index::

    git rm filename

* see diff between your branch (here master) and another one (here the remote origin)::

    git diff origin master

* see the list of files which differ::

    git diff origin master --stat

* propagate your changes to the main repository::

    git push


Binaries
--------

Old binaries generated for some old platforms may be download from here: https://gforge.inria.fr/frs/?group_id=9
