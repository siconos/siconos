
resultfile ="NE_BouncingBeam.dat"
#resultfile ="NE_BouncingBeam_D1MinusLinear.dat"


#plot resultfile u 1:2
resultfile1 ="NE_BouncingBeam_beam.dat"
#resultfile1 ="NE_BouncingBeam_D1MinusLinear_beam.dat"
set xrange [0.2:1.1]
set label 1 sprintf("Siconos NewtonEulerDS simulation.")  at first 0.5 , first 0.5
i=0
set term aqua 0
#plot resultfile1 u 3*i+1:3*i+3 w lp t 'beam1'
set term aqua 1
plot  resultfile u 1:9 w l t 'y'
set term aqua 2
plot  resultfile u 1:10 w l t 'ydot'
set term gif animate
set term aqua 3
plot  resultfile u 1:11 w l t 'lambda'
set term aqua 4
plot  resultfile u 1:12 w l t 'lambda1'
set term gif animate

outputfile="animatedbouncingbeam.gif"
outputfile="animatedbouncingbeam_D1MinusLinear.gif"
set label 1 sprintf("Siconos NewtonEulerDS simulation.")  at first -0.5 , first 0.5

set xrange [-1.05:1.05]
set yrange [-3:-0.9]
set size square
set output outputfile
n=1000   #n frames
#n=50
i=0
every=5

load "animonebeam.gp"
