.. _running_siconos:

Running a simulation
====================

Basics
------

At this step, you must have :

* installed properly siconos as explained in :ref:`siconos_install_guide`
* created a directory for your simulation, 
* written in this directory a 'driver' file, either in C++ or in python, as presented in :ref:`siconos_tutorials` or in one of the numerous examples of :ref:`siconos_examples`

Then, simply use siconos script on you driver (driver.cpp or driver.py)::

  siconos driver.cpp

Reminder: if *siconos_install_path* is not a standard path of your system, you may need to set some environment variables, mainly:

    * append *siconos_install_path*/bin to PATH

  
Plugins mechanism
-----------------

Some classes (mainly for dynamical systems and relations description), proposed a plugin (callback) mechanism, to allow user to set its own function to compute some
specific operators. This is detailed in :ref:`siconos_plugins`.

A 'plugin' is a dynamic library, generated from c++ source file, providing a set of functions that could be dynamically called by siconos objects.

The rules are :

* Siconos will consider any directory named XXXPluginXXX (XXX being whatever you want) or 'plugins' as a plugin directory
* for each plugin directory, siconos will create a library named XXXPluginXXX.ext (ext depends on your system) from all sources files in this directory

See for instance Mechanics/BouncingBall where external forces are defined with a plugin.

Extra source files
------------------

If needed, some user-defined source files can be taken into account. This may be useful to define some new classes, derived from standard siconos classes
(see for instance example Biology/StepSystem), to interface siconos with an other software and so on.
Anyway, to do this, just run::

  siconos --src_dir=path_to_extra_src YourDriver.cpp/.py

where *path_to_extra_src* is the directory where extra source files are saved.
Notice that by default, siconos will always check for additional sources in "src" dir of the current directory.
