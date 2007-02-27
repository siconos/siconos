#----------------------------------------------------------------------------------
# File enabling the comparison in gnuplot between :
#    - reference results from the SMASH simulation
#      of the diode bridge circuit with backward Euler
#      integrator and fixed 1 us time step (file "SMASHN0p25BE1us.dat").
#      The first 217 us are discarded since they correspond to the initialization
#      of the oscillating circuit.
#    - results from SICONOS  "DiodeBridgeBE1us.dat" (rename of "DiodeBridge.dat")
#      assumed to be obtained by a simulation of the DiodeBridge.cpp
#      sample with theta = 1 (Backward Euler integration) and a
#      time step of 1 us.
#
# Displayed results are :
#    - the voltage across the capacitor (or inductor)
#    - the current through the inductor
#    - the potential of each pole of the resistor (reference = diode R1 anode)
#    - the current through the resistor
#
# Slight differences may appear due to the smooth model of the diode
# used by SMASH. With a value of N (emission coefficient) of 0.25, the
# diode is a stiff one, with a threshold of around 0.25 V. Such a threshold
# yields small differences in the conduction state of the diode with respect
# to the ideal diode behavior currently (12/09/2005) modeled in the Siconos sample.
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
 
set grid xtics

set multiplot

set size winratio,winheight

set origin 0.0,winheight*3.0+heightoff
set ylabel "V" 1
plot\
  "SMASHN0p25BE1us.dat" u ($1-toffset):2 t "capacitor voltage , SMASH   BE 1us" w l,\
  "DiodeBridgeBE1us.dat" u 1:2 t           "capacitor voltage , SICONOS BE 1us" w l 
set origin 0.0,winheight*2.0+heightoff
set ylabel "A" 1
plot\
  "SMASHN0p25BE1us.dat" u ($1-toffset):3 t "inductor current , SMASH   BE 1us" w l,\
  "DiodeBridgeBE1us.dat" u 1:3 t           "inductor current , SICONOS BE 1us" w l
set origin 0.0,winheight+heightoff
set ylabel "V" 1
plot\
  "SMASHN0p25BE1us.dat" u ($1-toffset):4 t   "diode R1 cathod pot. , SMASH   BE 1us" w l,\
  "DiodeBridgeBE1us.dat" u 1:(-($5)) t       "diode R1 cathod pot. , SICONOS BE 1us" w l,\
  "SMASHN0p25BE1us.dat" u ($1-toffset):5 t   "diode F2 anode pot. ,  SMASH   BE 1us" w l,\
  "DiodeBridgeBE1us.dat" u 1:6 t             "diode F2 anode pot. ,  SICONOS BE 1us" w l

set bmargin 0
set format
set xtics axis
set xlabel "time in s"
set origin 0.0,0.0+heightoff
set ylabel "A" 1
plot\
  "SMASHN0p25BE1us.dat" u ($1-toffset):6 t "resistor current , SMASH   BE 1us" w l,\
  "DiodeBridgeBE1us.dat" u 1:($4 + $7) t   "resistor current , SICONOS BE 1us" w l
set nomultiplot

