
resultfile ="NE_1DS_1Knee_MLCP.dat"
resultfile_ref ="NE_1DS_1Knee_MLCP.ref"
#resultfile ="NE_1DS_1Knee_MLCP_D1MinusLinear.dat"
#resultfile_ref ="NE_1DS_1Knee_MLCP_D1MinusLinear.ref"

#plot resultfile u 1:2
resultfile1 ="NE_1DS_1Knee_MLCP_beam1.dat"
#resultfile1 ="NE_1DS_1Knee_MLCP_D1MinusLinear_beam1.dat"

set label 1 sprintf("Siconos NewtonEulerDS simulation.")  at first 0.5 , first 0.5
i=0
set term aqua 0
plot resultfile1 u 3*i+1:3*i+3 w lp t 'beam1'
set term aqua 1
plot  resultfile u 1:9 w l t 'norm2(y)'
set term aqua 2
plot  resultfile u 1:10 w l t 'norm2(ydot)'
set term aqua 3
set xrange [.0:0.5]
plot  resultfile u 1:2 w l t 'q[1]',\
resultfile_ref every ::1  u 1:2 w l t 'q[1]'

set term gif animate

outputfile="animatedonebeam.gif"
outputfile="animatedonebeam_D1MinusLinear.gif"
set label 1 sprintf("Siconos NewtonEulerDS simulation.")  at first -0.5 , first 0.5

set xrange [-1.05:1.05]
set yrange [-2:0.1]
set size square
set output outputfile
n=1000   #n frames
#n=50
i=0
every=5

load "animonebeam.gp"
