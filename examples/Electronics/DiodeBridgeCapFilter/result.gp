#----------------------------------------------------------------------------------
# File enabling the comparison in gnuplot between :
#    - reference results from the SMASH simulation
#      of the filtered diode bridge circuit with backward Euler
#      integrator and fixed 1 us time step (file "SMASHN0p25BE1us.dat").
#      The first 217 us are discarded since they correspond
#      to the initialization of the oscillating circuit.
#      
#    - results from SICONOS  "DiodeBridgeCapFilterBE1us.dat" 
#      (rename of "DiodeBridgeCapFilter.dat") assumed to be obtained
#      by a simulation of the DiodeBridgeCapFilter.cpp
#      sample with theta = 1 (Backward Euler integration) and a
#      time step of 1 us.
#
# Displayed results are :
#    - the voltage across the oscillator
#    - the potential of each pole of the resistor (reference = diode R1 anode)
#    - the voltage across the resistor
#    - the current through two diodes of the bridge
#
# - Slight differences may appear due to the smooth model of the diode
# used by SMASH. With a value of N (emission coefficient) of 0.25, the
# diode is a stiff one, with a threshold of around 0.25 V. Such a threshold
# yields small differences in the conduction state of the diode with respect
# to the ideal diode behavior currently (12/09/2005) modelled in the Siconos sample.
# - Larger differences appear on the potentials of the bridge's load (2nd graph). 
# They are due to the phases when all diodes are blocked : the load is floating with respect
# to the oscillator, and there is an infinity of solutions (bounded by the oscillator poles
# potentials). Thus, these differences are not meaningful.
#----------------------------------------------------------------------------------

#set  term X11 


toffset = 2.170E-04

basheight = 0.2
heightoff = 0.2
winratio = 1.0
winheight = basheight*winratio

set lmargin 8 
set bmargin 0
set tmargin 0

set size 1.0 , 1.0

set xzeroaxis
set format x ""
set xrange[2E-6:0.005]

set mxtics 
set grid xtics
set grid mxtics

set multiplot

set size winratio,winheight

set origin 0.0,winheight*3.0+heightoff
set ylabel "V" 
plot \
  "SMASHN0p25BE1us.dat" u ($1-toffset):2 t     "capacitor voltage , SMASH   BE 1us" w l,\
  "DiodeBridgeCapFilter.dat" u 1:2 t      "capacitor voltage , SICONOS BE 1us" w l
set origin 0.0,winheight*2.0+heightoff
plot\
  "SMASHN0p25BE1us.dat" u ($1-toffset):5 t          "diode R1 cathod pot. , SMASH   BE 1us" w l,\
  "DiodeBridgeCapFilter.dat" u 1:(- $5) t      "diode R1 cathod pot. , SICONOS BE 1us" w l,\
  "SMASHN0p25BE1us.dat" u ($1-toffset):3 t          "diode F2 anode pot. ,  SMASH   BE 1us" w l,\
  "DiodeBridgeCapFilter.dat" u 1:6 t           "diode F2 anode pot. ,  SICONOS BE 1us" w l
set origin 0.0,winheight+heightoff
plot\
  "SMASHN0p25BE1us.dat" u ($1-toffset):4 t              "resistor voltage , SMASH   BE 1us" w l,\
  "DiodeBridgeCapFilter.dat" u 1:(-($5 + $6)) t    "resistor voltage , SICONOS BE 1us" w l

set bmargin 0
set format
set xtics axis
set xlabel "time in s"
set origin 0.0,0.0+heightoff
set ylabel "A" 
plot\
  "SMASHN0p25BE1us.dat" u ($1-toffset):6 t     "diode F1 current , SMASH   BE 1us" w l,\
  "DiodeBridgeCapFilter.dat" u 1:7 t      "diode F1 current , SICONOS BE 1us" w l,\
  "SMASHN0p25BE1us.dat" u ($1-toffset):7 t     "diode R1 current , SMASH   BE 1us" w l,\
  "DiodeBridgeCapFilter.dat" u 1:4 t      "diode R1 current , SICONOS BE 1us" w l
set nomultiplot

