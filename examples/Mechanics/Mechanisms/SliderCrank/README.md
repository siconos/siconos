Slider Crank
============

This is an example of a simple multibody system shown on the figure below:

![image](slider_crank.*)

[![A video fo this example](https://img.youtube.com/vi/Wj0ZMcESw-Y/0.jpg)](https://youtu.be/Wj0ZMcESw-Y)

Model Definition
----------------

The model is defined in bodydef.py. The example is based on CAD files located in ./CAD and use OpenCascade. (version 0.16 or 0.18).

The option WITH\_CLEARANCE\_ON\_RODE can be set to 1 to add clearance between the rode 1 and 2.

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