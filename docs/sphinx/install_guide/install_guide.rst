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
* Check :ref:`siconos_dependencies`. Most of them are commonly or at least easily installed
  on many standard systems.
  
The quick way
-------------

If you do not want to bother with all installations details and only need a 'standard' siconos install (after ensuring the prerequisite described above):

* Start your favorite Python environment (e.g. `Micromaba <https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html>`_, `venv <https://docs.python.org/3/library/venv.html>`_ ...)

* Run

.. code-block:: bash
  
  cmake -S path_to_sources -B build-siconos
  cmake --build build-siconos -j N
  cmake --install build-siconos -j N
  ctest --test-dir build-siconos # optional
  

* N being the number of processes available on your system to run the compilation.
* path-to-sources: full path to Siconos sources
* build-siconos: the binary dir (where compilation, link, build take place). Choose any name/path. It must be different from path-to-sources.
  
If all went fine, you will get a full siconos installation as detailed in :ref:`siconos_whatsinstalled`
If not, step to :ref:`siconos_detailed_install`.

To check the install, try

.. code-block:: bash

   siconos --info
   
.. _siconos_detailed_install:
   
Detailed installation
---------------------

The first step of the installation process consists in running 'cmake'

.. code-block:: bash

   cmake -S path_to_sources -B build-siconos -DOPTION1_NAME=option1_value  ...

This command will explore your system, to generate an appropriate configuration, and some makefiles for Siconos, taking into account
some extra options, set as shown above. Many extra options exists to customize your build/install process of siconos.

The easiest way to handle Siconos options is to save them in a config file and to call cmake like this

.. code-block:: bash

   cmake -S path_to_sources -B build-siconos -DUSER_OPTIONS_FILE=option_file.cmake

Examples of options files are available in the directory `config_samples <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/tree/master/config_samples?ref_type=heads>`_ of Siconos source dir. To write your own file, just copy the file default.cmake somewhere and modify it according to your needs.

Choose any place/name for build-siconos, the only requirement is that it must be different from path_to_sources. This is a temporary directory that can be removed once the installation is done.

.. note::
   In place of the command-line cmake, you can also run::

      ccmake path_to_sources ...

   to open some dialog-interface to cmake configuration. 'cmake-gui' is also another option. For details check cmake documentation : https://cmake.org/runningcmake/ .
   

Once the cmake process is done, generated files stay in *build-siconos*, including a Makefile and a CMakeCache.txt. The latter contains all
the variables set during configuration. Do not forget to check the screen output of cmake to be sure that everything went fine.

Then you are ready to build siconos libraries and binaries::

  cmake --build build-siconos -j N

Or if you want to build a single target::

  cd build-siconos
  make target_name -j N

All available targets are obtained with::

  make help

Optionnaly (if WITH_TESTING is ON), you can run tests to check you build. See :ref:`siconos_run_tests`.

The last step is the installation of all required libraries, headers and so on in the right place::

  
  cmake --install build-siconos -j N

By default, everything will be installed 

- in your python env if it exists ($CONDA_PREFIX or $VIRTUAL_ENV)
- in $HOME/.siconos if not. In that case, add $HOME/.siconos/bin to your PATH so that siconos command can be found, e.g.::

    export PATH=$HOME/.siconos/bin:$PATH
  
Run

.. code-block:: bash

   siconos --info

to collect information about Siconos installation.

.. note::
   By default, no root privileged are required to run Siconos installation.

   We strongly recommend to retain the default installation mode, but there are other options:

   * use option SICONOS_INSTALL_SYSTEM_WIDE=true to install the software in the standard paths (/usr/local ...).
     This requires root privileges

   * use SICONOS_CUSTOM_INSTALL=<someplace> to customize Siconos installation path. Siconos binaries and libs will go to <someplace>.
     In that case, if ISOLATED_INSTALL=false (default) Siconos Python packages will remain in the default python install path.
     If ISOLATED_INSTALL=true, everything (including Siconos Python packages) will be installed in <someplace>

.. _siconos_whatsinstalled:

What will be installed?
-----------------------

.. note::

   Check the output of cmake to get info on where things will be installed for your current config
   

We denote *siconos_install_path* as the Siconos install root dir

* default:
  
  * *siconos_install_path* = $CONDA_PREFIX or $VIRTUAL_ENV if they exist
  * *siconos_install_path* = $HOME/.siconos if not 

* *siconos_install_path* = <someplace> if SICONOS_CUSTOM_INSTALL=<someplace> option was used with cmake
    
Then, the following files will be installed:

* *siconos_install_path*/lib/ with all shared libraries of the siconos components you asked for.
* *siconos_install_path*/include/siconos/ with all headers files needed by siconos
* *siconos_install_path*/share/siconos/ : extra files like cmake configuration, doc or anything that may be required at runtime
* *siconos_install_path*/bin/siconos : a script to run siconos simulation (see :ref:`siconos_runexample`).
* Python Siconos packages:

  * By default in the current default Python site-package ($CONDA_PREFIX, $VIRTUAL_ENV or user site), e.g.
    $HOME/siconosenv/lib/python3.10/site-packages or $HOME/.local/lib/python3.10/site-packages

  * if SICONOS_CUSTOM_INSTALL=<someplace> and ISOLATED_INSTALL=true in
    someplace/lib/python3.XY/site-packages (XY being your Python version)


.. warning::
   
   if *siconos_install_path* is not a standard path of your system, you may need to set some environment variables, mainly:

   * append *siconos_install_path*/bin to PATH
   * append path to Siconos Python packages to PYTHONPATH


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

  cd build-siconos
  make -j  test

To run only a set of tests, for example number 10 to 14::

  ctest -VV -I 10,14

'-V' or '-VV' is used to enable verbose and extra verbose mode. For other options, try 'man ctest' or check ctest documentation, https://cmake.org/documentation/.

To get a list of all available tests::

  ctest -N

To run python tests only::

  cd build-siconos
  py.test

Or in verbose mode::
  
  cd build-siconos
  py.test -s -v

Just a specific python test::
  
  cd build-siconos
  py.test -s -v wrap/siconos/tests/test_lcp.py

Concerning py.test, see http://pytest.org/latest/ or::
  py.test -h


Developers or advanced users options
------------------------------------
  

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
	
* WITH_BULLET=ON/OFF : enable/disable bullet (http://bulletphysics.org/wordpress/) for contact detection.

  Bullet minimal required version is 3.17.
  
  WITH_BULLET can be replaced by Bullet_ROOT=<some_path> to specify the path to your Bullet installation.

  Moreover, if you don't want to bother you with Bullet install, run

  .. code-block:: bash

	cmake -DBULLET_INSTALL=ON ...

  Bullet will be downloaded, built and installed as a siconos part, at the same place as Siconos.

  Last option, you can use the script ci_gitlab/Dockerfiles/install_bullet.sh to install Bullet 3.21 on your system (need to be root or sudo).

  .. code-block:: bash

        export CI_PROJECT_DIR=<some path where bullet will be cloned and built>
	source ci_gitlab/Dockerfiles/install_bullet.sh 


* WITH_OpenCASCADE=ON/OFF : enable/disable OpenCascade bindings (https://github.com/tpaviot/pythonocc-core)

.. note::
  for most of the required or optional dependencies, you can add some hints regarding their installation path to
  help cmake find them by using the option 'XXX_ROOT=<install_path>', XXX being the name of the package to be searched.
  For example::

    cmake -DFCLIB_ROOT=... 

.. _siconos_runexample:

Test your installation
----------------------

When all the installation process is done, you can test your installation by running a simple example.
(for non-standard installation path, mind :ref:`siconos_install_note`.). Try one of the numerous files
provided in `Siconos Tutorial project<https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/examples>`_


.. code-block:: bash

   git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials
   cd siconos-tutorials/examples/mechanics/BouncingBall

   siconos BouncingBallTS.cpp

   

You can also test all examples in a raw::

  git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials
  cmake -S siconos-tutorials/examples -B build-examples
  cmake --build build-examples -jN
  ctest --test-dir build-examples 

This will compile, link and execute all the examples distributed with siconos.

Check :ref:`running_siconos` for more details on *siconos* script.


Repositories
============
Install Siconos using the official repositories. We provide packages
for the distributions listed below.

.. warning::

   The packages below are not frequently updated and the Siconos versions available might be outdated.

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
  
