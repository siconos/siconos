.. _siconos_install_guide:

Installation guide
==================

Whatever your system is, you will need first to :

* Download the sources of Siconos as explained in :ref:`siconos_gettingsources`.
* Create anywhere (but not in *path_to_sources*) a build directory.
* Check :ref:`siconos_dependencies`. Most of them are commonly or at least easily installed
  on many standard systems.
* Check the list of :ref:`siconos_valid_platforms` on which siconos has been succesfully installed.
The quick way
-------------
If you do not want to bother with all installations details and only need a 'standard' siconos install, just do ::

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

    *to open some dialog-interface to cmake configuration. 'cmake-gui' is also another option. For details check cmake documentation : https://cmake.org/runningcmake/.*

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
 
.. _siconos_gettingsources:

Getting the sources
-------------------

Siconos is freely available on github, https://github.com/siconos/siconos, so for example just run::

   git clone git@github.com:siconos/siconos.git

and you will get a directory 'siconos' which contains all the sources. This directory will be denoted *path_to_sources* in the following.


.. _siconos_package:

Siconos package description
---------------------------
Siconos software is made of 4 components:

* numerics (C api). A collection of low-level algorithms for solving basic Algebra and optimization problem arising in the simulation of nonsmooth dynamical systems.

* kernel (C++ api), used to model and simulate nonsmooth dynamical systems.

* control (C++ api)

* mechanics (C++ api)


.. image:: /figures/siconos_components.*
	   
TODO : describe siconos distribution (main directories, files and so on)
  

.. _siconos_run_tests:

Running siconos tests
---------------------

You must enable tests with option WITH_TESTING=ON for cmake. To activate tests only for some chosen component, use::

  cmake -DWITH_<COMPONENT_NAME>_TESTING=ON

Then to run all tests::

  make -j N test

To run only a set of tests, for example number 10 to 14::

  ctest -VV -I 10,14

'-V' or '-VV' is used to enable verbose and extra verbose mode. For other options, try 'man ctest' or check ctest documentation, https://cmake.org/documentation/.

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
--------------------------------------

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

Developers or advanced users options
""""""""""""""""""""""""""""""""""""
  
* DEV_MODE=ON (OFF) : activate developper mode, which means for example some more aggressive options for compilations, more outputs and so on

* WITH_MUMPS=ON/OFF : to enable/disable mumps library (http://mumps.enseeiht.fr)

* WITH_FCLIB=ON/OFF : to enable/disable fclib interface

* WITH_DOXYGEN_WARNINGS=ON/OFF : verbose mode to explore doxygen warnings generated for siconos

* WITH_SERIALIZATION :

* WITH_GENERATION:

* WITH_CXX=ON/OFF : to enable/disable c++ compilation of the numerics package

* BUILD_SHARED_LIBS=ON/OFF : to build shared (ON) or static (OFF) for the siconos package.

* WITH_BULLET=ON/OFF : enable/disable bullet (http://bulletphysics.org/wordpress/) for contact detection.

* WITH_OCC=ON/OFF : enable/disable OpenCascade bindings (https://github.com/tpaviot/oce)

* WITH_FREECAD=ON/OFF : enable/disable Freecad python bindings (http://www.freecadweb.org)

* WITH_MECHANISMS=ON/OFF: enable/disable usage of Saladyn machanisms toolbox.

* WITH_DOXY2SWIG=ON/OFF : enable/disable conversion of doxygen outputs to python docstrings

For example, to build siconos with documentation for all components, no python bindings and an installation in '/home/myname/mysiconos', just run::

  cd build_directory
  cmake -DCMAKE_INSTALL_PREFIX='/home/myname/mysiconos' -DWITH_PYTHON_WRAPPER=OFF -DWITH_DOCUMENTATION=ON *path_to_sources*

But when you need a lot of options, this may get a bit tedious, with very long command line. To avoid this, you can use
:ref:`siconos_install_with_user_options`.

.. _siconos_install_with_user_options:

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

  # --------- User-defined options ---------
  # Use cmake -DOPTION_NAME=some-value ... to modify default value.
  # !!! Warning : do not suppress any line below, just set ON/OFF value !!!
  option(WITH_DOCUMENTATION "Build Documentation. Default = OFF" ON)
  option(WITH_PYTHON_WRAPPER "Build python bindings using swig. Default = ON" ON)
  option(WITH_DOXYGEN_WARNINGS "Explore doxygen warnings." OFF)
  option(WITH_DOXY2SWIG "Build swig docstrings from doxygen xml output. Default = ON." OFF)
  option(WITH_SYSTEM_INFO "Verbose mode to get some system/arch details. Default = off." OFF)
  option(WITH_TESTING "Enable 'make test' target" OFF)
  option(WITH_GIT "Consider sources are under GIT" OFF)
  option(WITH_SERIALIZATION "Compilation of serialization functions. Default = OFF" OFF)
  option(WITH_GENERATION "Generation of serialization functions with gccxml. Default = OFF" OFF)
  option(WITH_CXX "Enable CXX compiler for Numerics. Default=ON." ON)
  option(WITH_UNSTABLE "Enable this to include all 'unstable' sources. Default=OFF" OFF)
  option(BUILD_SHARED_LIBS "Building of shared libraries" ON)
  option(DEV_MODE "Compilation flags setup for developpers. Default: ON" OFF)
  option(WITH_BULLET "compilation with Bullet Bindings. Default = OFF" OFF)
  option(WITH_OCC "compilation with OpenCascade Bindings. Default = OFF" OFF)
  option(WITH_MUMPS "Compilation with MUMPS solver. Default = OFF" OFF)
  option(WITH_FCLIB "link with fclib when this mode is enable. Default = off." OFF)
  option(WITH_FREECAD "Use FreeCAD" OFF)
  option(WITH_MECHANISMS "Generation of bindings for Saladyn Mechanisms toolbox" OFF)
  option(WITH_XML "Enable xml files i/o. Default = ON" ON)
  # Set python install mode:
  # - user --> behave as 'python setup.py install --user'
  # - standard --> install in python site-package (ie behave as python setup.py install)
  # - prefix --> install in python CMAKE_INSTALL_PREFIX (ie behave as python setup.py install --prefix=CMAKE_INSTALL_PREFIX)
  set(siconos_python_install "user" CACHE STRING "Install mode for siconos python package")
  # List of components to build and installed
  # List of siconos component to be installed
  # complete list = externals numerics kernel control mechanics io
  set(COMPONENTS externals numerics kernel CACHE INTERNAL "List of siconos components to build and install")

  

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
