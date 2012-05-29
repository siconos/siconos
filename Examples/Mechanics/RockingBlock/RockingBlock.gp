set autoscale
set term X11
set term pdf

l=1;
L=1.5;


set xlabel 'time (s)' 
#set xrange [0.0001:10]   
#set yrange [0.0001:10]
set grid 

scheme =  "MoreauTSFixed"
#scheme  =  "ProjectedMoreauTS"
#scheme = "CombinedProjectedMoreauTS"

timestep= "150"

set output scheme."-RockingBlock-q-".timestep.".pdf"

plot 'result.dat'  using 1:2  t 'x' w lines,\
     'result.dat'  using 1:3  t 'y' w lines,\
     'result.dat'  using 1:4  t 'theta' w lines

set output scheme."-RockingBlock-v-".timestep.".pdf"

plot 'result.dat'  using 1:5  t 'dotx' w lines,\
     'result.dat'  using 1:6  t 'doty' w lines,\
     'result.dat'  using 1:7  t 'dottheta' w lines

set output scheme."RockingBlock-Energy-".timestep.".pdf"
plot 'result.dat'  using 1:8  t 'Kinetic' w lines,\
     'result.dat'  using 1:9  t 'Potential' w lines,\
     'result.dat'  using ($1):($8+$9)  t 'Total' w lines

set output scheme."-RockingBlock-g-".timestep.".pdf"

plot 'result.dat'  using 1:($3-1/2.0*(L*cos($4) -l*sin($4) ))  t 'g1' w lines,\
     'result.dat'  using 1:($3-1/2.0*(L*cos($4) +l*sin($4) ))  t 'g2' w lines

set output scheme."-RockingBlock-viol-".timestep.".pdf"

 plot 'result.dat'  using 1:  ( ($3-1/2.0*(L*cos($4) -l*sin($4) ) < 0 ) ?  ($3-1/2.0*(L*cos($4) -l*sin($4) ) )  : 0.0 )  t 'g1' w lines,\
      'result.dat'  using  1:  ( ($3-1/2.0*(L*cos($4) +l*sin($4) ) < 0 ) ?  ($3-1/2.0*(L*cos($4) +l*sin($4) ) )  : 0.0 ) t 'g2' w lines

