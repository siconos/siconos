Siconos mechanisms examples
---------------------------

[![A video of a circuit breaker simulation](https://img.youtube.com/vi/BvPUsuGX2jo/0.jpg)](https://youtu.be/BvPUsuGX2jo)


This directory contains some tests ot the mechanisms toolbox using  CAD files and OpenCascade.

Model Definition
----------------

The multibody system is defined in bodydef.py.

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
