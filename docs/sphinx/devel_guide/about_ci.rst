.. _about_ci:

Continuous integration for Siconos project
==========================================

A git push into siconos github will launch travis (on github) and jenkins (ci-inria) process, as described in .travis.yml file (for both travis and jenkins).

The following process will be executed on a worker:

* git clone siconos from github
* create build dir
* run in build dir::
    
    ../CI/driver.py --run --root-dir=..


For Travis, the 'worker' is the default one, provided by github-travis service.
The task executed is 'default' defined in tasks.py.
The ouput can be checked here : https://travis-ci.org/siconos/siconos

For Jenkins, the 'worker' is one of the 'nodes' available on Jenkins interface, e.g. siconos---vm0, siconos---vm1 ...
The executed tasks are those listed in tasks.py in the dictionnary known_tasks, for the chosen hostname (i.e. node name).
For example, if::

  known_tasks = {'siconos---vm0':
               (siconos_fedora_latest,
                siconos_gcc_asan,
     
               'siconos---vm1':
               (siconos_documentation,
                siconos_numerics_only,
     

then tasks siconos_fedora_latest and siconos_gcc_asan will be executed on node vm0, while siconos_documentation and siconos_numerics_only
will be executed on vm1.

Example of output: https://ci.inria.fr/siconos--/job/continuous3/lastSuccessfulBuild/console



Notice that this is equivalent to run, on another worker::

   ./driver.py  --run --root-dir=.. --tasks=siconos_documentation,siconos_numerics_only


To check what will be executed by the command above, just try::

   ./driver.py  --dry-run --run --root-dir=.. --tasks=siconos_documentation,siconos_numerics_only


The output will look like::

  Would call:
  - cmake -DMODE=Continuous -DCI_CONFIG=with_documentation -DWITH_DOCKER=1 -DBUILD_CONFIGURATION=Release -DDOCKER_DISTRIB=ubuntu:16.10 -DDOCKER_TEMPLATES=build-base,gcc,gfortran,gnu-c++,atlas-lapack,lpsolve,python-env,documentation -DDOCKER_TEMPLATE=gcc-atlas-lapack-documentation -DDOCKER_PROJECT_SOURCE_DIR=/home/perignon/Softs/siconos/. -DDOCKER_SHARED_DIRECTORIES= /home/perignon/Softs/siconos/./CI
  - make - ki target,
  for target in docker-build docker-cmake docker-make docker-make-install docker-make-doc docker-make-upload
  both from path _ubuntu-16.10_with_documentation

  Would call:
  - cmake -DMODE=Continuous -DCI_CONFIG=no_cxx -DWITH_DOCKER=1 -DBUILD_CONFIGURATION=Release -DDOCKER_DISTRIB=ubuntu:16.10 -DDOCKER_TEMPLATES=build-base,gcc,gfortran,atlas-lapack,lpsolve,python-env -DDOCKER_TEMPLATE=gcc-atlas-lapack -DDOCKER_PROJECT_SOURCE_DIR=/home/perignon/Softs/siconos/. -DDOCKER_SHARED_DIRECTORIES= /home/perignon/Softs/siconos/./CI
  - make - ki target,
  for target in docker-build docker-ctest
  both from path _ubuntu-16.10_no_cxx


Which means that (for numerics_only), in dir _ubuntu-16.10_no_cxx

* the 'cmake -DMODE ...' line will be executed to configure the ci project and its targets

The following targets will be executed:
* make docker-build, to create a docker image, based on ubuntu 16.10 with packages gcc, atlas ... will be created, with some associated volumes
* make docker-ctest : call cmake, make, make test on Siconos sources inside a docker container, based on the image created with make docker-build.

A report will be sent to siconos cdash.




For details on the description of a task, check :ref:`adding_ci_task`.
