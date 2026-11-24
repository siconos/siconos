.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _siconos_install_guide:

Build and install
#################

Overview
========

Once built and installed, the *Siconos software* consists of dynamic libraries, C++ headers, a few scripts and Python packages.
To use it, you need to fulfil some prerequisites, then download, build and install the package from the sources.

Binary packages for Debian are also available.



Installation from sources
=========================

Prerequisites
-------------

Whatever your system, you first need to:

* download the Siconos sources as explained in :ref:`download`;
* check the :ref:`siconos_dependencies`. Most of them are already installed, or at least easy to install,
  on standard systems.

The quick way
-------------

If you do not want to bother with the installation details and only need a 'standard' Siconos installation (once the prerequisites above are met):

* Start your favourite Python environment (e.g. `Micromamba <https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html>`_, `venv <https://docs.python.org/3/library/venv.html>`_ ...)

* Run

.. code-block:: bash

  cmake -S path_to_sources -B build-siconos
  cmake --build build-siconos -j N
  cmake --install build-siconos -j N
  ctest --test-dir build-siconos # optional


* N: the number of processes used to compile.
* path_to_sources: full path to the Siconos sources.
* build-siconos: the binary directory (where compilation and linking take place). Choose any name/path, as long as it differs from path_to_sources.

If all goes well, you get a full Siconos installation, as detailed in :ref:`siconos_whatsinstalled`.
If not, go to :ref:`siconos_detailed_install`.

To check the install, try

.. code-block:: bash

   siconos --info

.. _siconos_detailed_install:

Detailed installation
---------------------

The first step of the installation process is to run 'cmake':

.. code-block:: bash

   cmake -S path_to_sources -B build-siconos -DOPTION1_NAME=option1_value  ...

This command explores your system and generates a suitable configuration and makefiles for Siconos, taking into account
the extra options set as shown above. Many options exist to customize the build and installation of Siconos.

The easiest way to handle Siconos options is to save them in a configuration file and to call cmake like this:

.. code-block:: bash

   cmake -S path_to_sources -B build-siconos -DUSER_OPTIONS_FILE=option_file.cmake

Examples of options files are available in the directory `config_samples <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/tree/main/config_samples?ref_type=heads>`_ of the Siconos sources. To write your own file, copy the file default.cmake somewhere and modify it to suit your needs.

Choose any location/name for build-siconos; the only requirement is that it differs from path_to_sources. This is a temporary directory that can be removed once the installation is done.

.. note::
   Instead of the command-line cmake, you can also run::

      ccmake path_to_sources ...

   which opens an interactive interface for the cmake configuration. 'cmake-gui' is another option. For details, see the cmake documentation: https://cmake.org/runningcmake/.


Once cmake has run, the generated files are in *build-siconos*, including a Makefile and a CMakeCache.txt. The latter contains all
the variables set during the configuration. Check the cmake screen output to make sure that everything went fine.

You are then ready to build the Siconos libraries and binaries::

  cmake --build build-siconos -j N

Or if you want to build a single target::

  cd build-siconos
  make target_name -j N

All available targets are obtained with::

  make help

Optionally (if WITH_TESTING is ON), you can run the tests to check your build. See :ref:`siconos_run_tests`.

The last step installs all the required libraries, headers and so on in the right place::


  cmake --install build-siconos -j N

By default, everything is installed

- in your Python environment if there is one ($CONDA_PREFIX or $VIRTUAL_ENV);
- otherwise in $HOME/.siconos. In that case, add $HOME/.siconos/bin to your PATH so that the siconos command can be found, e.g.::

    export PATH=$HOME/.siconos/bin:$PATH

Run

.. code-block:: bash

   siconos --info

to get information about the Siconos installation.

.. note::
   By default, the Siconos installation does not require root privileges.

   We strongly recommend keeping the default installation mode, but there are other options:

   * use the option SICONOS_INSTALL_SYSTEM_WIDE=true to install the software in the standard paths (/usr/local...).
     This requires root privileges;

   * use SICONOS_CUSTOM_INSTALL=<someplace> to choose the Siconos installation path. The Siconos binaries and libraries go to <someplace>.
     In that case, if ISOLATED_INSTALL=false (default), the Siconos Python packages stay in the default Python installation path.
     If ISOLATED_INSTALL=true, everything (including the Siconos Python packages) is installed in <someplace>.

Selecting components
~~~~~~~~~~~~~~~~~~~~

Siconos is split into several components. By default, the standard components
(`externals`, `numerics`, `kernel`, `control`, `mechanics` and `io`) are built.

You can select components using ``WITH_<component>`` CMake options, for example::

  cmake -S path_to_sources -B build-siconos \
        -DWITH_numerics=ON -DWITH_kernel=ON -DWITH_mechanics=ON

The available components are:

* ``externals``
* ``numerics``
* ``kernel``
* ``control``
* ``mechanics``
* ``mechanisms``
* ``io``
* ``vkernel``

Dependencies are resolved automatically: enabling ``mechanics`` also enables
``kernel``, ``numerics`` and ``externals``.

Alternatively, you can still set the ``COMPONENTS`` list variable directly
(either on the command line or in a config file). When ``COMPONENTS`` is set,
it takes precedence over the ``WITH_<component>`` options.


.. _siconos_whatsinstalled:

What will be installed?
-----------------------

.. note::

   Check the cmake output to see where things will be installed with your current configuration.


We denote by *siconos_install_path* the Siconos installation root directory:

* default:

  * *siconos_install_path* = $CONDA_PREFIX or $VIRTUAL_ENV if they exist
  * *siconos_install_path* = $HOME/.siconos otherwise

* *siconos_install_path* = <someplace> if the SICONOS_CUSTOM_INSTALL=<someplace> option was given to cmake

The following files are then installed:

* *siconos_install_path*/lib/: all the shared libraries of the Siconos components you selected;
* *siconos_install_path*/include/siconos/: all the header files needed by Siconos;
* *siconos_install_path*/share/siconos/: extra files like the cmake configuration, docs, or anything that may be required at runtime;
* *siconos_install_path*/bin/siconos: a script to run Siconos simulations (see :ref:`siconos_runexample`);
* the Siconos Python packages:

  * by default in the current Python site-packages ($CONDA_PREFIX, $VIRTUAL_ENV or user site), e.g.
    $HOME/siconosenv/lib/python3.10/site-packages or $HOME/.local/lib/python3.10/site-packages;

  * if SICONOS_CUSTOM_INSTALL=<someplace> and ISOLATED_INSTALL=true, in
    someplace/lib/python3.XY/site-packages (XY being your Python version).


.. warning::

   If *siconos_install_path* is not a standard path on your system, you may need to set some environment variables, mainly:

   * append *siconos_install_path*/bin to PATH;
   * append the path to the Siconos Python packages to PYTHONPATH.


.. _siconos_package:

Siconos package description
---------------------------
The Siconos software is made of the components described below:

* **externals**: APIs or tools related to the external software libraries used by Siconos.

* **numerics** (C and Python API): a collection of low-level algorithms for solving the basic algebra and optimization problems arising in the simulation of nonsmooth dynamical systems.

* **kernel** (C++ and Python API): modeling and simulation of nonsmooth dynamical systems.

* **control** (C++ and Python API): control toolbox.

* **mechanics** (C++ and Python API): toolbox for collision detection and joints.

* **mechanisms** (C++ and Python API): toolbox for collision detection and joints (legacy version, will not be maintained in the long term).

* **io** (C++ API): tools related to input/output (HDF5, VTK...).


.. image:: /figures/siconos_components.*

The list of components to install can be set as shown in :ref:`siconos_detailed_install` (mind the dependencies shown in the figure above).


.. _siconos_run_tests:

Running the Siconos tests
-------------------------

To enable the tests, use the option WITH_TESTING=ON when running cmake.

Then, to run all the tests::

  cd build-siconos
  make -j  test

To run only a set of tests, for example numbers 10 to 14::

  ctest -VV -I 10,14

'-V' or '-VV' enables the verbose or extra verbose mode. For other options, try 'man ctest' or see the ctest documentation, https://cmake.org/documentation/.

To get the list of all available tests::

  ctest -N

To run the Python tests only::

  cd build-siconos
  py.test

Or in verbose mode::

  cd build-siconos
  py.test -s -v

A single Python test::

  cd build-siconos
  py.test -s -v wrap/siconos/tests/test_lcp.py

For py.test, see http://pytest.org/latest/ or::

  py.test -h


Options for developers and advanced users
-----------------------------------------


* WARNINGS_LEVEL: sets the level of compiler diagnostics.

  * 0: no warnings (default)
  * 1: enables many standard warnings (Wall, Wextra...). This should be the setup for developers.
  * 2: strict level: warnings become errors, and so on.

* WITH_MUMPS=ON/OFF: enables/disables the MUMPS library (http://mumps.enseeiht.fr).

* WITH_FCLIB=ON/OFF: enables/disables the fclib interface (https://github.com/FrictionalContactLibrary/fclib).
  This option is ON by default.
  The latest version of fclib (master branch of the GitHub repository) is downloaded and installed automatically as part of Siconos.
  If you need a specific version, or prefer to use a version already installed on your system, add the following option to your cmake command:

  .. code-block:: bash

	cmake -DFCLIB_ROOT=<path-to-your-fclib-installation> ...

* WITH_BULLET=ON/OFF: enables/disables Bullet (http://bulletphysics.org/wordpress/) for contact detection.

  The minimum required version of Bullet is 3.17.

  WITH_BULLET can be replaced by Bullet_ROOT=<some_path> to give the path to your Bullet installation.

  If you do not want to bother with installing Bullet, run

  .. code-block:: bash

	cmake -DBULLET_INSTALL=ON ...

  Bullet is then downloaded, built and installed as part of Siconos, in the same place as Siconos.

  Finally, you can use the script ci_gitlab/Dockerfiles/install_bullet.sh to install Bullet 3.21 on your system (requires root or sudo).

  .. code-block:: bash

        export CI_PROJECT_DIR=<some path where bullet will be cloned and built>
	source ci_gitlab/Dockerfiles/install_bullet.sh


* WITH_OpenCASCADE=ON/OFF: enables/disables the OpenCASCADE bindings (https://github.com/tpaviot/pythonocc-core).

.. note::
  For most of the required or optional dependencies, you can help cmake find them by giving their installation path
  with the option 'XXX_ROOT=<install_path>', XXX being the name of the package.
  For example::

    cmake -DFCLIB_ROOT=...

.. _siconos_runexample:

Test your installation
----------------------

Once the installation is done, you can test it by running a simple example.
Try one of the many files provided in the `Siconos tutorials project <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/>`_


.. code-block:: bash

   git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials
   cd siconos-tutorials/examples/mechanics/BouncingBall

   siconos BouncingBallTS.cpp



You can also test all the examples in a row::

  git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials
  cmake -S siconos-tutorials/examples -B build-examples
  cmake --build build-examples -jN
  ctest --test-dir build-examples

This compiles, links and runs all the examples distributed with Siconos.

See :ref:`running_siconos` for more details on the *siconos* script.


Repositories
============

Siconos can be installed from the official repositories of the distributions listed below.

.. warning::

   The packages below are not updated often, and the available Siconos versions may be outdated.

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
