set autoscale
#set term X11
set term pdf monochrome dashed

l=1;
L=1.5;


set xlabel 'time (s)' 
#set xrange [0.0001:10]   
#set yrange [0.0001:10]
set grid 

#scheme =  "Moreau"
scheme  = "Moreau-ProjectOnConstraints"
#scheme = "Moreau-CombinedProjection"
scheme = "D1MinusLinear"

timestep= "1e-5"
model= "-SliderCrank-Lagrangian-"

set output scheme.model.timestep.".pdf"
resultfile =  'result.dat'
#resultfile = 'Moreau-SliderCrank-Lagrangian-1e-6.dat'
#resultfile = 'Moreau-SliderCrank-Lagrangian-1e-5.dat'
#resultfile = 'Moreau-SliderCrank-Lagrangian-1e-4.dat'
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
plot resultfile  using ($2):($6) w lines  notitle

set auto
set xlabel "connecting-rod angle (rad)"
set ylabel "connecting rod speed(rad/s)"
plot resultfile  using  ($3):($6) w lines  notitle

# set xlabel "crank revolutions"
# set ylabel "X-Slider position"
# plot resultfile  using ($2):($12) w lines notitle
# set ylabel  "Y-Slider position"
# plot resultfile  using ($2):($12) w lines t notitle


set xlabel "X-Slider position"
set ylabel "Y-Slider position"

set xrange [-1.1:1.1]
set yrange [-1.1:1.1]

plot resultfile  using 12:13   w lines notitle

unset xrange
unset yrange
set auto
set xlabel 'Crank revolutions' 

basheight = 0.2
heightoff = 0.11
winratio = 1.0
winheight = basheight*winratio
hoffset = 0.05
set lmargin 8 
set bmargin 0
set tmargin 0

set size 1.0 , 1.0

set xzeroaxis
set format x ""
unset xlabel
#set grid xtics
set multiplot

set size winratio-0.1,winheight

set ytics -1.0,0.5,1.
set yrange [-1.25:1.25]
set xrange [0:2.0]
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin hoffset,(winheight+0.01)*3.0+heightoff
set ylabel "corner 1"
plot resultfile  using 2:8   notitle w lines
set origin hoffset,(winheight+0.01)*2.0+heightoff
set ylabel "corner 2"
plot resultfile  using 2:9  notitle w lines
set origin hoffset,(winheight+0.01)*1.0+heightoff
set ylabel "corner 3"
plot resultfile  using 2:10  notitle w lines
set bmargin 0
set format
set xtics border
set xlabel "Crank revolutions"
set ylabel "corner 4"
set origin hoffset,winheight*0.0+heightoff
plot resultfile  using 2:11  notitle w lines
     
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

set auto
set xrange [0:2.0]

set multiplot
set grid xtics

set ytics 0.0e-3,0.5e-3,2.0e-3
set yrange [-0.00025:0.00225]
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin hoffset,winheight*3.0+heightoff
plot resultfile  using 2:14  t 'g1' w lines
set origin hoffset,winheight*2.0+heightoff
plot resultfile  using 2:15  t 'g2' w lines
set origin hoffset,winheight*1.0+heightoff
plot resultfile  using 2:16  t 'g3' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin hoffset,winheight*0.0+heightoff
plot resultfile  using 2:17  t 'g4' w lines


unset multiplot
set xzeroaxis
set format x ""
unset xlabel
set grid xtics
set auto
set ytics axis
set xtics axis
set multiplot
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin hoffset,winheight*3.0+heightoff
plot resultfile  using 2:18  t 'dot g1' w lines
set origin hoffset,winheight*2.0+heightoff
plot resultfile  using 2:19  t 'dot g2' w lines
set origin hoffset,winheight*1.0+heightoff
plot resultfile  using 2:20  t 'dot g3' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin hoffset,winheight*0.0+heightoff
plot resultfile  using 2:21  t 'dot g4' w lines



unset multiplot
set xzeroaxis
set format x ""
unset xlabel
set grid xtics

set multiplot
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin hoffset,winheight*3.0+heightoff
plot resultfile  using 2:22  t 'lambda1' w lines
set origin hoffset,winheight*2.0+heightoff
plot resultfile  using 2:23  t 'lambda2' w lines
set origin hoffset,winheight*1.0+heightoff
plot resultfile  using 2:24  t 'lambda3' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin hoffset,winheight*0.0+heightoff
plot resultfile  using 2:25  t 'lambda4' w lines


unset multiplot
set xzeroaxis
set format x ""
unset xlabel
set grid xtics
set auto
set multiplot
set xrange [0:0.02]
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin hoffset,winheight*3.0+heightoff
plot resultfile  using 1:27  t 'active contact' w lines
set origin hoffset,winheight*2.0+heightoff
plot resultfile  using 1:22  t 'lambda1' w lines
set origin hoffset,winheight*1.0+heightoff
plot resultfile  using 1:23  t 'lambda2' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin hoffset,winheight*0.0+heightoff
plot resultfile  using 1:4  t 'slider angle' w lines



unset multiplot
reset

set auto

plot resultfile  using 1:4  t 'slider angle' w lines
set xrange [0:0.2]
set yrange [-0.5:2.5]
plot resultfile  using 1:27  t 'active contact' w lines
set auto
set xrange [0:0.25]
set ylabel "U_1"
set xlabel "Crank revolutions"
plot resultfile  using 2:18  notitle w lines
set xrange [0.15:0.25]
set ylabel "U_1"
set xlabel "Crank revolutions"
plot resultfile  using 2:18  notitle w lines





# set xzeroaxis
# set format x ""
# unset xlabel
# set grid xtics

# set multiplot
# #set output scheme."-SliderCrank-v-".timestep.".pdf"
# set origin hoffset,winheight*3.0+heightoff
# plot resultfile  using 2:26  t 'lambda(0) 1' w lines
# set origin hoffset,winheight*2.0+heightoff
# plot resultfile  using 2:27  t 'lambda(0) 2' w lines
# set origin hoffset,winheight*1.0+heightoff
# plot resultfile  using 2:28  t 'lambda(0) 3' w lines
# set bmargin 0
# set format
# set xtics axis
# set xlabel "Crank revolutions"
# set origin hoffset,winheight*0.0+heightoff
# plot resultfile  using 2:29  t 'lambda(0) 4' w lines



# reset
# set xzeroaxis
# set format x ""
# unset xlabel
# set grid xtics

# set multiplot
# #set output scheme."-SliderCrank-v-".timestep.".pdf"
# set origin hoffset,winheight*3.0+heightoff
# plot resultfile  using 2:14  t 'g1' w lines
# set origin hoffset,winheight*2.0+heightoff
# plot resultfile  using 2:18  t 'dot g1' w lines
# set origin hoffset,winheight*1.0+heightoff
# plot resultfile  using 2:22  t 'lambda1' w lines
# set bmargin 0
# set format
# set xtics axis
# set xlabel "Crank revolutions"
# set origin hoffset,winheight*0.0+heightoff
# plot resultfile  using 2:26  t 'lambda(0) 1' w lines


