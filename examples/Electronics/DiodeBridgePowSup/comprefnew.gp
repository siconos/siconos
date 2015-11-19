#----------------------------------------------------------------------------------
# File enabling the comparison in gnuplot between :
#    - reference results from SICONOS "DiodeBridgePowSupTRAP1us.dat".
#      
#    - recent results from SICONOS  "DiodeBridgePowSup.dat" 
#      obtained with a new version of SICONOS and the DiodeBridgePowSup.cpp
#      sample with theta = 0.5 (trapezoidal integration) and a
#      time step of 1 us.
#
# Displayed results are :
#    - the supply voltage 
#    - the current through two diodes of the bridge
#    - the voltage across the resistor
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
set xrange[3e-6:1.5e-3]

set mxtics 
set grid xtics
set grid mxtics

set multiplot

set size winratio,winheight

set origin 0.0,winheight*3.0+heightoff
set ylabel "V"
plot \
  "DiodeBridgePowSupTRAP1us.dat" u 1:2 t      "Supply voltage , SICONOS ref TRAP 1us" w l,\
  "DiodeBridgePowSup.dat"        u 1:2 t      "Supply voltage , SICONOS new TRAP 1us" w l

set origin 0.0,winheight*2.0+heightoff
set ylabel "A"
plot\
  "DiodeBridgePowSupTRAP1us.dat" u 1:7 t      "diode F1 current , SICONOS ref TRAP 1us" w l,\
  "DiodeBridgePowSup.dat"        u 1:7 t      "diode F1 current , SICONOS new TRAP 1us" w l

set origin 0.0,winheight+heightoff
plot\
  "DiodeBridgePowSupTRAP1us.dat" u 1:4 t      "diode R1 current , SICONOS ref TRAP 1us" w l,\
  "DiodeBridgePowSup.dat"        u 1:4 t      "diode R1 current , SICONOS new TRAP 1us" w l

set bmargin 0
set format
set xtics axis
set xlabel "time in s"
set origin 0.0,0.0+heightoff
set ylabel "V"
plot\
  "DiodeBridgePowSupTRAP1us.dat" u 1:(-($5 + $6)) t    "resistor voltage , SICONOS ref TRAP 1us" w l,\
  "DiodeBridgePowSup.dat"        u 1:(-($5 + $6)) t    "resistor voltage , SICONOS new TRAP 1us" w l

set nomultiplot

