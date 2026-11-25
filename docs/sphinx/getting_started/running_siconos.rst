.. SPDX-FileCopyrightText: INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _running_siconos:

Running a simulation
====================

Basics
------

At this point, you should have:

* installed Siconos as explained in :ref:`siconos_install_guide`,
* created a directory for your simulation,
* written a 'driver' file in this directory, in C++ or in Python, as in the many examples of :ref:`siconos_examples`.

Then run the siconos script on your driver (driver.cpp or driver.py)::

  siconos driver.cpp

Reminder: if *siconos_install_path* is not a standard path on your system, you may need to set some environment variables, mainly:

    * append *siconos_install_path*/bin to PATH


Plugins mechanism
-----------------

Some classes (mainly dynamical systems and relations) provide a plugin (callback) mechanism that lets users set their own functions to compute some
specific operators. This is detailed in :ref:`siconos_plugins`.

A 'plugin' is a dynamic library, built from C++ source files, providing a set of functions that Siconos objects can call dynamically.

The rules are:

* Siconos considers any directory named XXXPluginXXX (XXX being anything you want) or 'plugins' as a plugin directory;
* for each plugin directory, Siconos creates a library named XXXPluginXXX.ext (ext depends on your system) from all the source files in this directory.

See for instance Mechanics/BouncingBall, where external forces are defined with a plugin.

Extra source files
------------------

If needed, user-defined source files can be taken into account. This is useful to define new classes derived from standard Siconos classes
(see for instance the example Biology/StepSystem), to interface Siconos with other software, and so on.
To do this, run::

  siconos --src_dir=path_to_extra_src YourDriver.cpp/.py

where *path_to_extra_src* is the directory containing the extra source files.
Note that by default, Siconos always looks for additional sources in the "src" directory of the current directory.
