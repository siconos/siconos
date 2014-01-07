#----------------------------------------------------------------------------------
# File enabling the comparison in gnuplot between :
#    - reference results from the SMASH simulation
#      of the CircuitRLCD circuit with backward Euler
#      integrator and fixed 1 us time step (file "SMASHN0p25BE1us.dat").
#      The first 217 us are discarded since they correspond to the initialization
#      of the oscillating circuit.
#    - results from SICONOS  "CircuitRLCDBE1us.dat" (rename of "CircuitRLCD.dat")
#      assumed to be obtained by a simulation of the CircuitRLCD.cpp
#      sample with theta = 1 (Backward Euler integration) and a
#      time step of 1 us.
#
# Displayed results are :
#    - the voltage across the capacitor (or inductor)
#    - the current through the inductor
#    - the voltage across the resistor
#    - the current through the resistor
#
# Slight differences may appear due to the smooth model of the diode
# used by SMASH. With a value of N (emission coefficient) of 0.25, the
# diode is a stiff one, with a threshold of around 0.25 V. Such a threshold
# yields small differences in the conduction state of the diode with respect
# to the ideal diode behavior currently (12/09/2005) modeled in the Siconos sample.
#----------------------------------------------------------------------------------

#set  term X11 

#!tail -n 10000 CircuitRLCD.dat > resultgp.dat

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






set term post
set output "StateInputOutputFR.ps"
set multiplot

set size winratio,winheight

set origin 0.0,winheight*3.0+heightoff
set ylabel "V"
plot\
  "CircuitRLCD-gammaFR0.dat" u ($1-toffset):2 t "capacitor voltage , Gamma0   BE 1us" w l,\
  "CircuitRLCD-gammaFR1.dat" u ($1-toffset):2 t "capacitor voltage , Gamma1   BE 1us" w l,\
  "CircuitRLCD.dat" u 1:2 t      "capacitor voltage , Gamma1/2 BE 1us" with lines 
set origin 0.0,winheight*2.0+heightoff
set ylabel "A" 
plot\
  "CircuitRLCD-gammaFR0.dat" u ($1-toffset):3 t "inductor current , Gamma0   BE 1us" w l,\
  "CircuitRLCD-gammaFR1.dat" u ($1-toffset):3 t "inductor current , Gamma1   BE 1us" w l,\
  "CircuitRLCD.dat" u 1:3 t      "inductor current , Gamma1/2 BE 1us" w l
set origin 0.0,winheight+heightoff
set ylabel "V"
plot\
  "CircuitRLCD-gammaFR0.dat" u ($1-toffset):4 t "diode voltage , Gamma0   BE 1us" w l,\
  "CircuitRLCD-gammaFR1.dat" u ($1-toffset):4 t "diode voltage , Gamma1   BE 1us" w l,\
  "CircuitRLCD.dat" u 1:4 t      "diode voltage , Gamma1/2 BE 1us" w l,\
  "CircuitRLCD.dat" u 1:( -(4*pi*pi)  *$2 - 10*$5) t      "diode voltage(comp) , Gamma1/2 BE 1us" w l

set bmargin 0
set format
set xtics axis
set xlabel "time in s"
set origin 0.0,0.0+heightoff
set ylabel "A" 
plot\
  "CircuitRLCD-gammaFR0.dat" u ($1-toffset):5 t "diode current , Gamma0   BE 1us" w l,\
  "CircuitRLCD-gammaFR1.dat" u ($1-toffset):5 t "diode current , Gamma1   BE 1us" w l,\
  "CircuitRLCD.dat" u 1:5 t      "diode current , Gamma1/2 BE 1us" w l
set nomultiplot

set output "EnergyFR.ps"
set size
set origin

plot\
  "CircuitRLCD-gammaFR0.dat" u ($1-toffset):7 t "storage function , Gamma0   BE 1us" w l,\
  "CircuitRLCD-gammaFR1.dat" u ($1-toffset):7 t "storage function , Gamma1   BE 1us" w l,\
  "CircuitRLCD.dat" u 1:7 t      "storage function , Gamma1/2 BE 1us" w l

 
plot\
  "CircuitRLCD-gammaFR0.dat" u ($1-toffset):8 t "dissipation function , Gamma0   BE 1us" w l,\
  "CircuitRLCD-gammaFR1.dat" u ($1-toffset):8 t "dissipation function , Gamma1   BE 1us" w l,\
  "CircuitRLCD.dat" u 1:8 t "dissipation function , Gamma1/2 BE 1us" w l

plot\
  "CircuitRLCD-gammaFR0.dat" u ($1-toffset):9 t "total function , Gamma0   BE 1us" w l,\
  "CircuitRLCD-gammaFR1.dat" u ($1-toffset):9 t "total function , Gamma1   BE 1us" w l,\
  "CircuitRLCD.dat" u 1:9    t  "total function , Gamma1/2 BE 1us" w l


