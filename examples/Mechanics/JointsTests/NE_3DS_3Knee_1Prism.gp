
resultfile ="NE_3DS_3Knee_1Prism_GMP.dat"
#resultfile ="NE_3DS_3Knee_1Prism_MLCP_D1MinusLinear.dat"

resultfile ="NE_3DS_3Knee_1Prism_MLCP_MoreauJeanCombinedProjection.dat"
set xrange [-0.1:3]
set yrange [-3:0.1]
#plot resultfile u 1:2
resultfile1 ="NE_3DS_3Knee_1Prism_beam1.dat"
resultfile2 ="NE_3DS_3Knee_1Prism_beam2.dat"
resultfile3 ="NE_3DS_3Knee_1Prism_beam3.dat"
#resultfile1 ="NE_3DS_3Knee_1Prism_MLCP_D1MinusLinear_beam1.dat"
#resultfile2 ="NE_3DS_3Knee_1Prism_MLCP_D1MinusLinear_beam2.dat"
#resultfile3 ="NE_3DS_3Knee_1Prism_MLCP_D1MinusLinear_beam3.dat"
set label 1 sprintf("Siconos NewtonEulerDS simulation.")  at first 1.5 , first 0
i=0
plot resultfile1 u 3*i+1:3*i+3 w lp t 'beam1',  resultfile2 u 3*i+1:3*i+3 w lp t 'beam2',  resultfile3 u 3*i+1:3*i+3 w lp t 'beam3' 
set term gif animate

outputfile="animatedbeams.gif"
#outputfile="animatedbeams_D1MinusLinear.gif"

set size square
set xrange [-0.1:3]
set yrange [-3:0.1]
set output outputfile
n=1000    #n frames
n=500
i=0
every=5

load "animbeam.gp"
