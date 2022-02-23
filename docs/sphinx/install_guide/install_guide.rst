.. _siconos_install_guide:

Build and install
#################

.. contents::
   :local:


Overview
========

*Siconos software*, once properly build and installed, consists in some dynamics libraries, c++ headers, a few scripts and python packages.
To use the software, you need to fulfill some prerequesites, download build and install the package from sources.

Binaries for Debian are also available.

 

Installation from sources
=========================
      
Prerequisites
-------------

Whatever your system is, you will need first to :

* Download the sources of Siconos as explained in :ref:`download`.
* Create anywhere (but not in *path_to_sources*) a build directory.
* Check :ref:`siconos_dependencies`. Most of them are commonly or at least easily installed
  on many standard systems.
  
The quick way
-------------

If you do not want to bother with all installations details and only need a 'standard' siconos install, just do (first ensure you follow the steps described above : download, check dependencies ...)::

  cd path_to_build 
  cmake path_to_sources
  make -j N
  make install


N being the number of processes available on your system to run the compilation. Note that you'll probably need to run the install
command as root.

If all went fine, you will get a full siconos installation in /usr/local, as detailed in :ref:`siconos_whatsinstalled`
If not, step to :ref:`siconos_detailed_install`.

.. _siconos_detailed_install:
   
Detailed installation
---------------------

The first step of the installation process consists in running 'cmake'::

   cd path_to_build
   cmake path_to_sources -DOPTION1_NAME=option1_value -DOPTION2_NAME=option2_value ...

This command will explore your system, to generate an appropriate configuration, and some makefiles, for siconos, taking into account
some extra options, set as shown above. Many extra options exists to customize your build/install process of siconos and most of the available options
are detailled in :ref:`siconos_cmake_options`.


    *Remark : In place of the command-line cmake, you can also run*::

      ccmake path_to_sources ...

    *to open some dialog-interface to cmake configuration. 'cmake-gui' is also another option. For details check cmake documentation : https://cmake.org/runningcmake/ .*

Once the cmake process is done, you will get many generated files in *path_to_build*, including a Makefile and a CMakeCache.txt. The latter contains all
the variables set during configuration. Do not forget to check the screen output of cmake to be sure that everything went fine.

Then you are ready to build siconos libraries and binaries::

  make -j N

Or if you want to build a single target::

  make target_name -j N

All available targets are obtained with::

  make help

Optionnaly (if WITH_TESTING is ON), you can run tests to check you build. See :ref:`siconos_run_tests`.

The last step is the installation of all required libraries, headers and so on in the right place::

  make install -j N

Use CMAKE_INSTALL_PREFIX option to choose the path for your installation. If not set a default path is chosen, usually /usr/local (that depends on your system).
 
.. _siconos_package:

Siconos package description
---------------------------
Siconos software is made of different components described below

* **externals** : API or tools related to external software libraries used by Siconos.

* **numerics** (C  and Python api). A collection of low-level algorithms for solving basic Algebra and optimization problem arising in the simulation of nonsmooth dynamical systems.

* **kernel** (C++ and Python api), used to model and simulate nonsmooth dynamical systems.

* **control** (C++ and Python api) : control toolbox

* **mechanics** (C++ and Python api) : toolbox for collision detection and joints

* **mechanisms** (C++ and Python  api) : toolbox for collision detection and joints (legacy version, won't be sustained in long term)

* **io** (C++ api) : tools related to input/outputs (hdf5, vtk ...)


.. image:: /figures/siconos_components.*

The list of components to be installed can be set using :ref:`siconos_install_with_user_options` (mind the dependencies shown in the figure above).


.. _siconos_run_tests:

Running siconos tests
---------------------

To enable tests, use the option WITH_TESTING=ON when running cmake.

Then to run all tests::

  make -j  test

To run only a set of tests, for example number 10 to 14::

  ctest -VV -I 10,14

'-V' or '-VV' is used to enable verbose and extra verbose mode. For other options, try 'man ctest' or check ctest documentation, https://cmake.org/documentation/.

To get a list of all available tests::

  ctest -N

To run python tests only::

  cd path_to_build
  py.test

Or in verbose mode::
  
  cd path_to_build
  py.test -s -v

Just a specific python test::
  
  cd path_to_build
  py.test -s -v wrap/siconos/tests/test_lcp.py

Concerning py.test, see http://pytest.org/latest/ or::
  py.test -h

  
.. _siconos_whatsinstalled:

What will be installed?
-----------------------

For *siconos_install_path* being the value you choose for siconos install, running 'make install' will result in:


* *siconos_install_path*/lib/ with all shared libraries of the siconos components you asked for.
* *siconos_install_path*/include/siconos/ with all headers files needed by siconos
* *siconos_install_path*/share/siconos/ : extra files like cmake configuration, doc or anything that may be required at runtime
* *siconos_install_path*/bin/siconos : a script to run siconos simulation (see :ref:`siconos_runexample`).

.. _siconos_install_note:

Remark
""""""
if *siconos_install_path* is not a standard path of your system, you may need to set some environment variables, mainly:

* append *siconos_install_path*/bin to PATH


.. _siconos_cmake_options:

CMake options
-------------

Most options are like '-DWITH_XXX=ON or OFF to enable or disable some behavior or some interface to other libraries.
If ON, the cmake system will search for XXX libraries, headers, or anything required on your system and will end up in error if not found. 

Most common options
"""""""""""""""""""

* CMAKE_INSTALL_PREFIX=some_path : to change the default path of Siconos installation. Default depends on your system. For example on unix-like
  system, it is usually /usr/local.

* WITH_DOCUMENTATION=ON (OFF) : to enable (disable) the generation of siconos source code documentation and manuals generation.

* WITH_PYTHON_WRAPPER=ON (OFF) : to enable (disable) the generation of a python interface to siconos.

* WITH_CMAKE_BUILD_TYPE=Debug, Release, ... : to choose the build mode, i.e. the default compiler flags used to build siconos.

* WITH_TESTING : to enable/disable tests


.. warning::
  Sometimes you may try to modify some options without success.
  This because, when activated during the first run of cmake, some options are persistent. The only way to modify them is to remove the content of
  the build directory.

Developers or advanced users options
""""""""""""""""""""""""""""""""""""
  

* WARNINGS_LEVEL: to set compiler diagnostics level.

  * 0: no warnings (default)
  * 1: activate many standard warnings (Wall, Wextras ...). This should be the setup for developers.
  * 2: strict level, turn warnings to errors and so on.

* WITH_MUMPS=ON/OFF : to enable/disable mumps library (http://mumps.enseeiht.fr)

* WITH_FCLIB=ON/OFF : to enable/disable fclib interface (https://github.com/FrictionalContactLibrary/fclib). 
  This option is ON by default.
  The last version of fclib (master branch of the github repository) will be downloaded and installed automatically as part of Siconos
  If you need a specific version or prefer using a version already installed on your system, add the following option to your cmake command:

  .. code-block:: bash

	cmake -DFCLIB_ROOT=<path-to-your-fclib-installation> ...
	
* WITH_DOXYGEN_WARNINGS=ON/OFF : verbose mode to explore doxygen warnings generated for siconos

* WITH_SERIALIZATION:

* WITH_GENERATION:

* WITH_CXX=ON/OFF : to enable/disable c++ compilation of the numerics package. Must be ON if kernel component is used.

* BUILD_SHARED_LIBS=ON/OFF : to build shared (ON) or static (OFF) for the siconos package.

* WITH_BULLET=ON/OFF : enable/disable bullet (http://bulletphysics.org/wordpress/) for contact detection.

* WITH_OCE=ON/OFF : enable/disable OpenCascade bindings (https://github.com/tpaviot/oce)

* WITH_FREECAD=ON/OFF : enable/disable Freecad python bindings (http://www.freecadweb.org)

* WITH_DOXY2SWIG=ON/OFF : enable/disable conversion of doxygen outputs to python docstrings

For example, to build siconos with documentation for all components, no python bindings and an installation in '/home/myname/mysiconos', just run

.. code-block:: bash

  cd build_directory
  cmake -DCMAKE_INSTALL_PREFIX='/home/myname/mysiconos' -DWITH_PYTHON_WRAPPER=OFF -DWITH_DOCUMENTATION=ON *path_to_sources*

But when you need a lot of options, this may get a bit tedious, with very long command line. To avoid this, you can use :ref:`siconos_install_with_user_options`.

.. _siconos_install_with_user_options:

.. note::
  for most of the required or optional dependencies, you can add some hints regarding their installation path to
  help cmake find them by using the option 'XXX_ROOT=<install_path>', XXX being the name of the package to be searched.
  For example::

    cmake -DFCLIB_ROOT=... 

User-defined option file
------------------------

To avoid very long and boring command line during cmake call, you can write a 'myoption.cmake' and call::

  cd build_directory
  cmake -DUSER_OPTIONS_FILE=myoption.cmake path_to_sources

Warnings:

* your file MUST have the '.cmake' extension
* if you provide only its name to USER_OPTIONS_FILE, your file must be either in *path_to_sources* or in *path_to_build* directory
  else, you must give the absolute path to your file, for example::
     
    cmake -DUSER_OPTIONS_FILE=/home/myname/myoptions_for_siconos.cmake path_to_sources

To write your own file, just copy the file default_options.cmake (in *path_to_sources*/cmake) and modify it according to your needs.

Here is an example, to build numerics and kernel, with documentation, no tests ...::

  # --- List of siconos components to build and install ---
  # The complete list is : externals numerics kernel control mechanics mechanisms io
  set(COMPONENTS externals numerics kernel CACHE INTERNAL "List of siconos components to build and install")

  option(WITH_PYTHON_WRAPPER "Build and install python bindings using swig. Default = ON" ON)
  option(WITH_SERIALIZATION "Compilation of serialization functions. Default = OFF" OFF)
  option(WITH_GENERATION "Generation of serialization functions with doxygen XML. Default = OFF" OFF)

  # --- Build/compiling options ---
  set(WARNINGS_LEVEL 0 CACHE INTERNAL "Set compiler diagnostics level. 0: no warnings, 1: developer's minimal warnings, 2: strict level, warnings to errors and so on. Default =0")
  option(WITH_CXX "Enable CXX compiler for numerics. Default = ON" ON)
  option(WITH_FORTRAN "Enable Fortran compiler. Default = ON" ON)
  option(FORCE_SKIP_RPATH "Do not build shared libraries with rpath. Useful only for packaging. Default = OFF" OFF)
  option(NO_RUNTIME_BUILD_DEP "Do not check for runtime dependencies. Useful only for packaging. Default = OFF" OFF)
  option(WITH_UNSTABLE_TEST "Enable this to include all 'unstable' test. Default=OFF" OFF)
  option(BUILD_SHARED_LIBS "Building of shared libraries. Default = ON" ON)
  option(WITH_SYSTEM_INFO "Verbose mode to get some system/arch details. Default = OFF." OFF)
  option(WITH_TESTING "Enable 'make test' target" OFF)
  option(WITH_GIT "If true, try to get info (commit sha ...) from siconos sources git repository." OFF)

  # --- Documentation setup ---
  option(WITH_DOCUMENTATION "Build Documentation. Default = OFF" ON)
  option(WITH_DOXYGEN_WARNINGS "Explore doxygen warnings. Default = OFF" OFF)
  option(WITH_DOXY2SWIG "Build swig docstrings from doxygen xml output. Default = OFF." ON)

  # --- List of external libraries/dependencies to be searched (or not) ---
  option(WITH_BULLET "compilation with Bullet Bindings. Default = OFF" OFF)
  option(WITH_OCE "compilation with OpenCascade Bindings. Default = OFF" OFF)
  option(WITH_MUMPS "Compilation with the MUMPS solver. Default = OFF" OFF)
  option(WITH_UMFPACK "Compilation with the UMFPACK solver. Default = OFF" OFF)
  option(WITH_SUPERLU "Compilation with the SuperLU solver. Default = OFF" OFF)
  option(WITH_SUPERLU_MT "Compilation with the SuperLU solver, multithreaded version. Default = OFF" OFF)
  option(WITH_FCLIB "link with fclib when this mode is enable. Default = ON" ON)
  option(WITH_FREECAD "Use FreeCAD. Default = OFF" OFF)
  option(WITH_RENDERER "Install OCE renderer. Default = OFF" OFF)
  option(WITH_SYSTEM_SUITESPARSE "Use SuiteSparse installed on the system instead of built-in CXSparse library. Default = ON" ON)
  option(WITH_XML "Enable xml files i/o. Default = OFF" OFF)

  # -- Installation setup ---
  # Set python install mode:
  # - user --> behave as 'python setup.py install --user'
  # - standard --> install in python site-package (ie behave as python setup.py install)
  # - prefix --> install in python CMAKE_INSTALL_PREFIX (ie behave as python setup.py install --prefix=CMAKE_INSTALL_PREFIX)
  set(siconos_python_install "standard" CACHE STRING "Install mode for siconos python package")

  # If OFF, headers from libraries in externals will not be installed.
  option(INSTALL_EXTERNAL_HEADERS "Whether or not headers for external libraries should be installed. Default=OFF" OFF)

  # If ON, internal headers will not be installed.
  option(INSTALL_INTERNAL_HEADERS "Whether or not headers for internal definitions should be installed. Default=OFF" OFF)

  
.. _siconos_runexample:

Test your installation
----------------------

When all the installation process is done, you can test your installation by running a simple example.
(for non-standard installation path, mind :ref:`siconos_install_note`.). Try one of the numerous files
provided in Siconos Examples package::

  siconos BouncingBallTS.cpp


You can also test all examples in a raw::

  cd another_build_directory
  cmake path_to_sources/Examples
  make -jN
  make test


This will compile, link and execute all the examples distributed with siconos.

Check :ref:`running_siconos` for more details on *siconos* script.


Repositories
============
Install Siconos using the official repositories. We provide packages
for the distributions listed below.


Debian bullseye 
----------------

(thanks to Steven Sinclair's work)

.. code-block:: bash

   apt install siconos

FreeBSD 
--------

(thanks to yurivict, yuri@freebsd) 

.. code-block:: bash

  pkg install siconos
  
