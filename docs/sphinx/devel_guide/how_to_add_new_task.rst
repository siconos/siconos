.. _adding_ci_task:

How to add and test a new continuous integration task for Siconos project
=========================================================================


What is a "Continuous integration task"?
----------------------------------------

* the definition of a context (operating system, libraries, a specific compiler ...)
* a specific siconos configuration (list of components, with or without python and so on)
* a list of 'targets' to be executed (like build siconos, run tests ...)

Once defined (see below) this task will be executed on a runner (Travis or Jenkins/ci-inria) everytime somebody push into siconos github repository.

Create a task
-------------

1. Create a file my_task_name.cmake in CI/config dir. This file will be used to set siconos configuration options and components to be built.
  For example, if you need to compile siconos externals and numerics components, without python wrappers in DEV_MODE, task_name.cmake should contain::

    set_option(DEV_MODE ON)
    set_option(WITH_PYTHON_WRAPPER OFF)
    set_components(externals;numerics)

  By default, variables not set in my_task_name.cmake are read from cmake/siconos_default.cmake.

2. Create a new entry in file task.py, e.g.::

   my_new_task = CiTask(
     ci_config='task_name',
     distrib='ubuntu:16.10',
     pkgs=['build-base', 'gcc', 'gfortran', 'gnu-c++', 'atlas-lapack', 'python-minimal'],
     srcs=['.'],
     targets={'.': ['docker-build', 'docker-ctest']})

 * ci_config is the name (without extension) of the file in directory CI/config where siconos configuration is described.
    Note that the file task_name.cmake must exist (and may be empty).
 * distrib is the name of a docker image with its version (from docker hub)
 * pkgs is the list of dependencies that must be installed in the docker image
 * srcs is the path to what must be configured/built (e.g. siconos/ or siconos/examples sources)
 * targets is the list of what will be executed. See the list of existing targets below (:ref:`ci_targets`).  

3. add the task into the dictionnary known_tasks in file tasks.py. The key in this dictionnary will defined the worker on which
   the task will be executed.

Options sent by driver
----------------------


* For Docker (i.e. used when cmake is called to build docker container, images, volumes ...)

  * WITH_DOCKER : sounds rather useless, always 1.
  * DOCKER_DISTRIB : name of the docker image used as source. Set with parameter 'distrib' of tasks in tasks.py
  * DOCKER_TEMPLATES : list of dependencies (e.g. lapack, gcc, ...)
  * DOCKER_TEMPLATE : ? The difference with DOCKER_TEMPLATES is not clear to me
  * DOCKER_PROJECT_SOURCE_DIR : path to Siconos sources.
  * DOCKER_SHARED_DIRECTORIES :
  * 

* For ctest

  * BUILD_CONFIGURATION : build type used by ctest (ie CMAKE_BUILD_TYPE). Release, Debug, Profiling (default).
    Set with parameter 'build_configuration' of tasks in tasks.py.

    *Note FP: used only when 'make docker-ctest' target is called while 'make docker-cmake' will keep default conf from siconos. This should be fixed?*
  * MODE :
  * CI_CONFIG
  


.. _ci_targets:

Available ci targets
--------------------
All possible targets are described/defined in CI/cmake/Docker.cmake file. SOFT corresponds to the input srcs of task class, e.g. path to siconos
CMakeLists.txt or to siconos/examples CMakeLists.txt.

* docker-build : create docker image and associated volumes (workdir and /usr/local)
* docker-cmake : configure SOFT inside docker image
* docker-make : build SOFT inside docker image
* docker-make-clean : clean build dir
* docker-make-install : install SOFT inside docker image
* docker-make-uninstall : uninstall SOFT
* docker-make-test : run SOFT tests (if any)
* docker-make-doc : build and publish siconos doc
* docker-ctest : run ctest for SOFT (i.e. cmake, make, make test)
* docker-hard-clean : clean docker on worker (remove unused volumes, images ...). Execute script docker-cleanup.sh.
* docker-interactive : start a docker container based on the created image
