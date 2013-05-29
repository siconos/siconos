
resultfile ="NE_3DS_3Knee_1Prism_MLCP.dat"
set xrange [-0.1:3]
set yrange [-3:0.1]
#plot resultfile u 1:2
resultfile1 ="NE_3DS_3Knee_1Prism_MLCP_beam1.dat"
resultfile2 ="NE_3DS_3Knee_1Prism_MLCP_beam2.dat"
resultfile3 ="NE_3DS_3Knee_1Prism_MLCP_beam3.dat"

plot resultfile1 u 1:3 w lp,  resultfile2 u 1:3 w lp,  resultfile3 u 1:3 w lp 

set term gif animate 
set output "animatedbeams.gif"
n=1000    #n frames
i=0
every=5

load "animbeam.gp"
