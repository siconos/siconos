#----------------------------------------------------------------------------------
# File enabling the comparison in gnuplot between :
#    - reference results from the SPICE simulation
#      of the filtered diode bridge circuit with trapezoidal
#      integrator and fixed 0.1 us time step (file "SPICEN0p25TRAP0p1us.dat").
#      
#    - results from SICONOS  "DiodeBridgePowSupTRAP0p1us.dat" 
#      obtained by a simulation of the DiodeBridgePowSup.cpp
#      sample with theta = 0.5 (trapezoidal integration) and a
#      time step of 0.1 us.
#
# Displayed results are :
#    - the supply voltage 
#    - the current through two diodes of the bridge
#    - the voltage across the resistor
#
# The threshold of the diodes was taken into account in the SICONOS sample,
# yielding very small differences in the conduction state of the diodes.
#
#----------------------------------------------------------------------------------

toffset = 0.0

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

#      The first 0.3 us are discarded due to a high current linked
#      to the initial loading of the tank capacitor.

set xrange[0.3e-6:1.5e-3]

set mxtics 
set grid xtics
set grid mxtics

set multiplot

set size winratio,winheight

set origin 0.0,winheight*3.0+heightoff
set ylabel "V"
plot \
  "DiodeBridgePowSupTRAP0p1us.dat" u 1:2 t      "Supply voltage , SICONOS TRAP 0p1us" w l,\
  "SPICEN0p25TRAP0p1us.dat" u ($1-toffset):2 t     "Supply voltage , SPICE   TRAP 0p1us" w l

set origin 0.0,winheight*2.0+heightoff
set ylabel "A"
plot\
  "DiodeBridgePowSupTRAP0p1us.dat" u 1:7 t      "diode F1 current , SICONOS TRAP 0p1us" w l,\
  "SPICEN0p25TRAP0p1us.dat" u ($1-toffset):6 t     "diode F1 current , SPICE   TRAP 0p1us" w l

set origin 0.0,winheight+heightoff
plot\
  "DiodeBridgePowSupTRAP0p1us.dat" u 1:4 t      "diode R1 current , SICONOS TRAP 0p1us" w l,\
  "SPICEN0p25TRAP0p1us.dat" u ($1-toffset):7 t     "diode R1 current , SPICE   TRAP 0p1us" w l


set bmargin 0
set format
set xtics axis
set xlabel "time in s"
set origin 0.0,0.0+heightoff
set ylabel "V"
plot\
  "DiodeBridgePowSupTRAP0p1us.dat" u 1:(-($5 + $6)) t    "resistor voltage , SICONOS TRAP 0p1us" w l,\
  "SPICEN0p25TRAP0p1us.dat" u ($1-toffset):4 t              "resistor voltage , SPICE   TRAP 0p1us" w l

set nomultiplot

