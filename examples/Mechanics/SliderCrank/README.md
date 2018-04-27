Slider Crank
============

This is an example of a simple multibody system, developed in [4] and shown in the figure below:

![image](slider_crank.*)

All variables names and parameters values are those from the paper cited above.

In siconos, this example is simulated in several different ways :

-  SliderCrankMoreauJeanOSI.cpp: simulation with the standard Moreau--Jean one-step integrator
-  SliderCrankD1MinusLinear.cpp: simulation with the Time discontinuous Galerkin integrator D1 minus written at the acceleration level [2]
-  SliderCrankD1MinusLinearVelocityLevel.cpp: simulation with the Time discontinuous Galerkin integrator D1 minus written at the velocity level [3]
-  SliderCrankMoreauJeanCombinedProjectionOSI.cpp: simulation with the Moreau--Jean one-step integrator with the combined projection algorithm for the projection at the position [1]
-  SliderCrankMoreauJeanDirectProjectionOSI.cpp: simulation with the Moreau--Jean one-step integrator with the direct projection algorithm for the projection at the position [1]

[1] Vincent Acary, Projected event-capturing time-stepping schemes for nonsmooth mechanical systems with unilateral contact and Coulomb’s friction, Computer Methods in Applied Mechanics and Engineering, Volume 256, 2013, Pages 224-250. https://doi.org/10.1016/j.cma.2012.12.012.

[2] Thorsten Schindler, Vincent Acary, Timestepping schemes for nonsmooth dynamics based on discontinuous Galerkin methods: Definition and outlook, Mathematics and Computers in Simulation,Volume 95,2014, Pages 180-199, https://doi.org/10.1016/j.matcom.2012.04.012.

[3] Schindler, T., Rezaei, S., Kursawe, J., Acary, V.
Half-explicit timestepping schemes on velocity level based on time-discontinuous Galerkin methods. Computer Methods in Applied Mechanics and Engineering, Volume 290, Issue undefined, June 05, 2015

[4] Flores, P., Leine, R. & Glocker, C. Multibody Syst Dyn (2010) 23: 165. https://doi.org/10.1007/s11044-009-9178-y

Usage
-----

    siconos example_name.cpp

and to plot the results, you can use for example:

    gnuplot -p result.gp

This must lead to Fig 11 (d) of the paper




