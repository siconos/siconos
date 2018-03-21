.. _slider_crank_example:

Slider Crank
============

.. highlight:: c++

This is an example of a simple multibody system, developed in :cite:`FloresGlockerLeine` and shown on the figure below:

.. image:: /figures/mechanics/slider_crank/slider_crank.*


In siconos, this example is simulated in several different ways :

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

or, for python files using mechanisms toolbox::

  cd examples/Mechanics/Mechanisms/
  siconos -P siconos-mechanism.py .

and to plot the results, you can use for example::

  gnuplot -p result.gp

This must lead to Fig 11 (d) of the paper

The example is based on cad files located in examples/Mechanics/Mechanisms/SliderCrank/CAD

The option WITH_CLEARANCE_ON_RODE can be set to 1 to add clearance between the rode 1 and 2.

All variables names and parameters values are those from the paper cited above.
