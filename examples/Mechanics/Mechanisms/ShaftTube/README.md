Slider Crank
============

This is an example of a simple multibody system shown in the video below:

[![A video fo this example](https://img.youtube.com/vi/hsqYWCo0Fu4/0.jpg)](https://youtu.be/hsqYWCo0Fu4W)

Model Definition
----------------

The model is defined in bodydef.py. The example is based on CAD files located in ./CAD and use OpenCascade. (version 0.16 or 0.18).

The option WITH\_CLEARANCE\_ON\_RODE can be set to 1 to add clearance between the rods 1 and 2.

The general options of the siconos_mechanisms toolbox can be overwritten in the file `mbtbLocalOptions.py`


GUI and visualization
---------------------

A light-weight GUI has been developped base on pythonocc (version 0.16 or 0.18).	

The option `with3D=1` in the file `mbtbLocalOptions.py` enables the 3D GUI.

Usage
-----

The simulation may be executed with this command:

    siconos -Psiconos_mechanisms .

or directly

    siconos_mechanisms

The simulation results can also be viewed from the hdf5 file :

    siconos_vview siconos_mechanisms.hdf5

