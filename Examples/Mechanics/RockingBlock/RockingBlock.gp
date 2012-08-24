set autoscale
set term X11
#set term pdf
#extension= ".pdf"
#set term tikz standalone color solid size 5in,3in
set term tikz standalone monochrome  size 5in,3in
extension=".tex"



l=1;
L=1.5;


set xlabel 'time (s)' 
#set xrange [0.0001:10]   
#set yrange [0.0001:10]
set grid 

#scheme =  "MoreauTSFixed"
#scheme  =  "ProjectedMoreauTS"
scheme = "CombinedProjectedMoreauTS"

timestep= "150"

set output scheme."-RockingBlock-q-".timestep.extension

plot 'result.dat'  using 1:2  t '$x$' w lines,\
     'result.dat'  using 1:3  t '$y$' w lines,\
     'result.dat'  using 1:4  t '$\theta$' w lines

set output scheme."-RockingBlock-v-".timestep.extension

plot 'result.dat'  using 1:5  t '$\dot x$' w lines,\
     'result.dat'  using 1:6  t '$\dot y$' w lines,\
     'result.dat'  using 1:7  t '$\dot \theta$' w lines

set output scheme."-RockingBlock-Energy-".timestep.extension
plot 'result.dat'  using 1:8  t 'Kinetic energy' w lines,\
     'result.dat'  using 1:9  t 'Potential energy' w lines,\
     'result.dat'  using ($1):($8+$9)  t 'Total energy' w lines

set output scheme."-RockingBlock-g-".timestep.extension

plot 'result.dat'  using 1:($3-1/2.0*(L*cos($4) -l*sin($4) ))  t 'g1' w lines,\
     'result.dat'  using 1:($3-1/2.0*(L*cos($4) +l*sin($4) ))  t 'g2' w lines

set output scheme."-RockingBlock-viol-".timestep.extension

 plot 'result.dat'  using 1:  ( ($3-1/2.0*(L*cos($4) -l*sin($4) ) < 0 ) ?  ($3-1/2.0*(L*cos($4) -l*sin($4) ) )  : 0.0 )  t 'g1' w lines,\
      'result.dat'  using  1:  ( ($3-1/2.0*(L*cos($4) +l*sin($4) ) < 0 ) ?  ($3-1/2.0*(L*cos($4) +l*sin($4) ) )  : 0.0 ) t 'g2' w lines

