set autoscale
set term X11
set term pdf

l=1;
L=1.5;


set xlabel 'time (s)' 
#set xrange [0.0001:10]   
#set yrange [0.0001:10]
set grid 

#scheme =  "Moreau"
scheme  = "Moreau-ProjectOnConstraints"
#scheme = "Moreau-CombinedProjection"

timestep= "1e-6"

set output scheme."-SliderCrank-".timestep.".pdf"
resultfile =  'result.dat'

#set output scheme."-SliderCrank-q-".timestep.".pdf"

plot resultfile  using 1:2  t 'crank revolution' w lines,\
     resultfile  using 1:3  t 'rod angle' w lines,\
     resultfile  using 1:4  t 'slide angle' w lines

set xrange[0:2]
set xlabel "crank revolutions"
set ylabel "crank speed(rad/s)"
plot resultfile  using ($2):($5) w lines notitle

set xlabel "crank revolutions"
set ylabel "connecting-rod speed(rad/s)"
plot resultfile  using ($2):($6) w lines t ' '

set auto
set xlabel "connecting-rod angle (rad)"
set ylabel "connecting rod speed(rad/s)"
plot resultfile  using  ($3):($6) w lines t ' '

set xlabel "crank revolutions"
set ylabel "X-Slider position"
plot resultfile  using ($2):($12) w lines notitle
set ylabel  "Y-Slider position"
plot resultfile  using ($2):($12) w lines t ' '


set xlabel "X-Slider position"
set ylabel "Y-Slider position"







set xrange [-1.1:1.1]
set yrange [-1.1:1.1]

plot resultfile  using 12:13  t 'phase portrait' w lines

unset xrange
unset yrange
set auto
set xlabel 'Crank revolutions' 

basheight = 0.2
heightoff = 0.1
winratio = 1.0
winheight = basheight*winratio

set lmargin 8 
set bmargin 0
set tmargin 0

set size 1.0 , 1.0

set xzeroaxis
set format x ""
 unset xlabel
set grid xtics

set multiplot

set size winratio-0.1,winheight

#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin 0.0,winheight*3.0+heightoff
plot resultfile  using 2:8  t 'Corner 1' w lines
set origin 0.0,winheight*2.0+heightoff
plot resultfile  using 2:9  t 'Corner 2' w lines
set origin 0.0,winheight*1.0+heightoff
plot resultfile  using 2:10  t 'Corner 3' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin 0.0,winheight*0.0+heightoff
plot resultfile  using 2:11  t 'Corner 4' w lines
     
unset multiplot
#set output scheme."SliderCrank-Energy-".timestep.".pdf"
# plot resultfile  using 1:8  t 'Kinetic' w lines,\
#      'result.dat'  using 1:9  t 'Potential' w lines,\
#      'result.dat'  using ($1):($8+$9)  t 'Total' w lines
#set output scheme."-SliderCrank-g-".timestep.".pdf"
set xzeroaxis
set format x ""
 unset xlabel
set grid xtics

set multiplot
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin 0.0,winheight*3.0+heightoff
plot resultfile  using 2:14  t 'g1' w lines
set origin 0.0,winheight*2.0+heightoff
plot resultfile  using 2:15  t 'g2' w lines
set origin 0.0,winheight*1.0+heightoff
plot resultfile  using 2:16  t 'g3' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin 0.0,winheight*0.0+heightoff
plot resultfile  using 2:17  t 'g4' w lines


unset multiplot
set xzeroaxis
set format x ""
 unset xlabel
set grid xtics

set multiplot
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin 0.0,winheight*3.0+heightoff
plot resultfile  using 1:14  t 'g1' w lines
set origin 0.0,winheight*2.0+heightoff
plot resultfile  using 1:15  t 'g2' w lines
set origin 0.0,winheight*1.0+heightoff
plot resultfile  using 1:16  t 'g3' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Time(s)"
set origin 0.0,winheight*0.0+heightoff
plot resultfile  using 1:17  t 'g4' w lines

unset multiplot
set xzeroaxis
set format x ""
unset xlabel
set grid xtics

set multiplot
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin 0.0,winheight*3.0+heightoff
plot resultfile  using 2:18  t 'dot g1' w lines
set origin 0.0,winheight*2.0+heightoff
plot resultfile  using 2:19  t 'dot g2' w lines
set origin 0.0,winheight*1.0+heightoff
plot resultfile  using 2:20  t 'dot g3' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin 0.0,winheight*0.0+heightoff
plot resultfile  using 2:21  t 'dot g4' w lines



unset multiplot
set xzeroaxis
set format x ""
unset xlabel
set grid xtics

set multiplot
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin 0.0,winheight*3.0+heightoff
plot resultfile  using 2:22  t 'lambda1' w lines
set origin 0.0,winheight*2.0+heightoff
plot resultfile  using 2:23  t 'lambda2' w lines
set origin 0.0,winheight*1.0+heightoff
plot resultfile  using 2:24  t 'lambda3' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin 0.0,winheight*0.0+heightoff
plot resultfile  using 2:25  t 'lambda4' w lines

unset multiplot
set xzeroaxis
set format x ""
unset xlabel
set grid xtics

set multiplot
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin 0.0,winheight*3.0+heightoff
plot resultfile  using 2:26  t 'lambda(0) 1' w lines
set origin 0.0,winheight*2.0+heightoff
plot resultfile  using 2:27  t 'lambda(0) 2' w lines
set origin 0.0,winheight*1.0+heightoff
plot resultfile  using 2:28  t 'lambda(0) 3' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin 0.0,winheight*0.0+heightoff
plot resultfile  using 2:29  t 'lambda(0) 4' w lines



unset multiplot
set xzeroaxis
set format x ""
unset xlabel
set grid xtics

set multiplot
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin 0.0,winheight*3.0+heightoff
plot resultfile  using 2:14  t 'g1' w lines
set origin 0.0,winheight*2.0+heightoff
plot resultfile  using 2:18  t 'dot g1' w lines
set origin 0.0,winheight*1.0+heightoff
plot resultfile  using 2:22  t 'lambda1' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin 0.0,winheight*0.0+heightoff
plot resultfile  using 2:26  t 'lambda(0) 1' w lines




