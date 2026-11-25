.. SPDX-FileCopyrightText: INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _about_ci:

Continuous integration for the Siconos project
==============================================

A git push to the Siconos GitHub repository launches the Travis (GitHub) and Jenkins (ci-inria) processes, as described in the .travis.yml file (for both Travis and Jenkins).

The following process is executed on a worker:

* git clone siconos from GitHub
* create a build directory
* run in the build directory::

    ../CI/driver.py --run --root-dir=..


For Travis, the 'worker' is the default one, provided by the GitHub-Travis service.
The task executed is 'default', defined in tasks.py.
The output can be checked here: https://travis-ci.org/siconos/siconos

For Jenkins, the 'worker' is one of the 'nodes' available in the Jenkins interface, e.g. siconos---vm0, siconos---vm1, etc.
The tasks executed are those listed in tasks.py in the dictionary known_tasks, for the chosen hostname (i.e. node name).
For example, if

.. code-block:: python

  known_tasks = {'siconos---vm0':
               (siconos_fedora_latest,
                siconos_gcc_asan,

               'siconos---vm1':
               (siconos_documentation,
                siconos_numerics_only,


then the tasks siconos_fedora_latest and siconos_gcc_asan are executed on node vm0, while siconos_documentation and siconos_numerics_only
are executed on vm1.

Example of output: https://ci.inria.fr/siconos--/job/continuous3/lastSuccessfulBuild/console



Note that this is equivalent to running, on another worker::

   ./driver.py  --run --root-dir=.. --tasks=siconos_documentation,siconos_numerics_only


To check what the command above will execute, try::

   ./driver.py  --dry-run --run --root-dir=.. --tasks=siconos_documentation,siconos_numerics_only


The output will look like::

  Would call:
  - cmake -DMODE=Continuous -DCI_CONFIG=with_documentation -DWITH_DOCKER=1 -DBUILD_CONFIGURATION=Release -DDOCKER_DISTRIB=ubuntu:16.10 -DDOCKER_TEMPLATES=build-base,gcc,gfortran,gnu-c++,lpsolve,python-env,documentation -DDOCKER_TEMPLATE=gcc-openblas-lapacke-documentation -DDOCKER_PROJECT_SOURCE_DIR=/home/perignon/Softs/siconos/. -DDOCKER_SHARED_DIRECTORIES= /home/perignon/Softs/siconos/./CI
  - make - ki target,
  for target in docker-build docker-cmake docker-make docker-make-install docker-make-doc docker-make-upload
  both from path _ubuntu-16.10_with_documentation

  Would call:
  - cmake -DMODE=Continuous -DCI_CONFIG=no_cxx -DWITH_DOCKER=1 -DBUILD_CONFIGURATION=Release -DDOCKER_DISTRIB=ubuntu:16.10 -DDOCKER_TEMPLATES=build-base,gcc,gfortran,lpsolve,python-env -DDOCKER_TEMPLATE=openblas-lapacke -DDOCKER_PROJECT_SOURCE_DIR=/home/perignon/Softs/siconos/. -DDOCKER_SHARED_DIRECTORIES= /home/perignon/Softs/siconos/./CI
  - make - ki target,
  for target in docker-build docker-ctest
  both from path _ubuntu-16.10_no_cxx


This means that (for numerics_only), in the directory _ubuntu-16.10_no_cxx,

* the 'cmake -DMODE ...' line is executed to configure the CI project and its targets.

Then the following targets are executed:

* make docker-build: creates a Docker image based on Ubuntu 16.10 with the packages gcc, etc., and some associated volumes;
* make docker-ctest: runs cmake, make and make test on the Siconos sources inside a Docker container based on the image created by make docker-build.

A report is sent to the Siconos CDash.
