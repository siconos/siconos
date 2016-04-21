#set  term postscript
#set output "result.ps"
#set term X11

toffset = 2.170E-04

basheight = 0.33
heightoff = 0.05
winratio = 0.95
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

set origin 0.0,winheight*2.0+heightoff
set ylabel 'm'


resultfile="BallNewtonEuler.dat"
resulfile ="BallNewtonEulerOnMovingPlane.dat"
plot \
resultfile u 1:2 t "Ball position" w l,\
resultfile u 1:8 t "Plane position " w l

set origin 0.0,winheight*1.0+heightoff
set ylabel "m/s" 
plot \
resultfile u 1:3 t "Ball Velocity" w l,\
resultfile u 1:9 t "Plane Velocity" w l

set bmargin 0
set format
set xtics axis
set xlabel "time in s"
set origin 0.0,0.0+heightoff
set ylabel "A" 

set ylabel "N s " 
plot \
resultfile u 1:4 t "contact impulse" w l,\
resultfile u 1:12 t "reaction impulse on plane" w l


set nomultiplot

