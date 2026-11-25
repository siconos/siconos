.. SPDX-FileCopyrightText: INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _slider_crank_example:

Slider Crank
============

.. highlight:: c++

This is an example of a simple multibody system, introduced in :cite:`FloresGlockerLeine` and shown in the figure below:

.. image:: /figures/mechanics/slider_crank/slider_crank.*


In Siconos, this example is simulated in several different ways:

* examples/Mechanics/SliderCrank/SliderCrankMoreauJeanOSI.cpp
* examples/Mechanics/SliderCrank/SliderCrankD1MinusLinear.cpp
* examples/Mechanics/SliderCrank/SliderCrankD1MinusLinearVelocityLevel.cpp
* examples/Mechanics/SliderCrank/SliderCrankMoreauJeanCombinedProjectionOSI.cpp
* examples/Mechanics/SliderCrank/SliderCrankMoreauJeanDirectProjectionOSI.cpp
* examples/Mechanics/Mechanisms/SliderCrank/bodyref.py

Usage
-----

::

   siconos example_name.cpp

or, for Python files using the mechanisms toolbox::

  cd examples/Mechanics/Mechanisms/
  siconos -P siconos-mechanism.py .

To plot the results, you can use for example::

  gnuplot -p result.gp

This should reproduce Fig. 11 (d) of the paper.

The example is based on CAD files located in examples/Mechanics/Mechanisms/SliderCrank/CAD.

The option WITH_CLEARANCE_ON_RODE can be set to 1 to add clearance between rods 1 and 2.

All variable names and parameter values are those of the paper cited above.
