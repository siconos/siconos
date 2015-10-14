set autoscale
#set term X11
#set term pdf

l1 = 0.1530
l2 = 0.3060
a = 0.05
b = 0.025
c = 0.001

set xlabel 'time (s)' 
#set xrange [0.0001:10]   
#set yrange [0.0001:10]
set grid 


# Define two helper functions
ismin(x) = (x<min)?min=x:0
ismax(x) = (x>max)?max=x:0

# Initialise the 'global' vars
max=-1e38
min=1e38

# Run through the data and pass it to the helper functions.
# Any expression for the 'using' will do, as long as it contains both
# helper functions
set output
#set term X11
plot "./simu.txt" u 1:(ismin($30)*ismax($30))

print max
print min
magnify = 1000



#scheme =  "Moreau"
#scheme  = "Moreau-ProjectOnConstraints"
scheme = "Moreau-CombinedProjection"
#model = "-SliderCrank-"
model = "-SliderCrank-r0e-0"
#model = "-SliderCrank-r6e-3"
timestep= "-h1e-4"
h= 1.e-4
set terminal pdf monochrome dashed

set output scheme.model.timestep.".pdf"
resultfile =  './simu.txt'
#resultfile =  'Moreau-SliderCrank-r0e-0-h1e-6.txt'
#set output scheme."-SliderCrank-q-".timestep.".pdf"

set xrange [0:0.25]
crankrevo(x,y,t) = (y<=0)?((t<0.1)?(acos(x)/pi):(acos(x)/pi+2.0)):(acos(-x)/pi+1.0)
crankangle(x,y,t) = (y<=0)?((t<0.1)?(2.0*acos(x)):(2.0*acos(x)+4.0*pi)):(2.0*acos(-x)+2.0*pi)
connectingangle(x,y,t) = (y<=0)?((2.0*acos(x))):(-2.0*acos(x))
sliderangle(x,y,t) = (y<=0)?((2.0*acos(x))):(-2.0*acos(x))

plot resultfile  using (h*$1):(crankrevo($5,$7,h*$1))  t 'crank revolution' w lines,\
     resultfile  using (h*$1):(connectingangle($19,$21,h*$21)) t 'rod angle' w lines,\
     resultfile  using (h*$1):(2.0*acos($33)) t 'slide angle' w lines
# plot resultfile  using (h*$1):(crankangle($5,$7,h*$1))  t 'crank angle' w lines
# plot resultfile  using (h*$1):(sliderangle($33,$35,h*$1))  t 'slider angle' w lines
#plot resultfile  using (h*$1):($5)  t 'q0' w lines,\
      resultfile  using (h*$1):($7)  t 'q2' w lines
#plot resultfile  using (h*$1):($33)  t 'q0' w lines,\
     resultfile  using (h*$1):($35)  t 'q2' w lines
#plot resultfile  using (h*$1):(($7<0)?(2.0*acos($5)/pi):(2.0*acos(-$5)/pi+2.0))  t 'crank revolution' w lines
#plot resultfile  using (h*$1):(crankrevo($5,$7))  t 'crank revolution' w lines

set xrange[0:2]
set xlabel "crank revolutions"
set ylabel "crank speed(rad/s)"
plot resultfile  using (crankrevo($5,$7,h*$1)):(-$13) w lines notitle


set xlabel "crank revolutions"
set ylabel "connecting-rod speed(rad/s)"
plot resultfile  using (crankrevo($5,$7,h*$1)):(-$27) w lines notitle

set auto
set xlabel "connecting-rod angle (rad)"
set ylabel "connecting rod speed(rad/s)"
plot resultfile  using  (connectingangle($19,$21,h*$21)):(-$27) w lines notitle

# set xlabel "crank revolutions"
# set ylabel "X-Slider position"
# plot resultfile  using (crankrevo($5,$7,h*$1)):($30) w lines notitle
# set ylabel  "Y-Slider position"
# plot resultfile  using (crankrevo($5,$7,h*$1)):($32) w lines t ' '


set xlabel "X-Slider position"
set ylabel "Y-Slider position"
set xrange[-1.05:0.05]
set yrange[-1.1:1.1]


plot resultfile  u (1/(max-min)*$30-1/(max-min)*min-1):(magnify*$32) w l notitle

set auto
set xlabel "time(s)"
set ylabel "Energy"
plot resultfile using (h*$1):($15)  t 'Kinetic energy crank' w lines,\
     resultfile using (h*$1):($29)  t 'Kinetic energy rod' w lines,\
     resultfile using (h*$1):($43)  t 'Kinetic energy slider' w lines,\
     resultfile using (h*$1):($15+$29+$43)  t 'Total Kinetic energy' w lines


unset xrange
unset yrange
set auto
set xrange[0:2.0]
set xlabel 'Crank revolutions' 
set ylabel ""
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

set size winratio,winheight

set ytics -1.0,0.5,1.
set yrange [-1.25:1.25]
#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin 0.0,(winheight+0.01)*3.0+heightoff
set ylabel "corner 1"
plot resultfile  using (crankrevo($5,$7,h*$1)):($32/c )  notitle w lines
set ylabel "corner 2"
set origin 0.0,(winheight+0.01)*2.0+heightoff
plot resultfile  using (crankrevo($5,$7,h*$1)):($32/c)  notitle w lines
set ylabel "corner 3"
set origin 0.0,(winheight+0.01)*1.0+heightoff
plot resultfile  using (crankrevo($5,$7,h*$1)):($32 /c) notitle w lines
set bmargin 0
set format
#set xtics axis
set ylabel "corner 4"
set xlabel "Crank revolutions"
set origin 0.0,winheight*0.0+heightoff
plot resultfile  using (crankrevo($5,$7,h*$1)):($32   /c)  notitle w lines


set xzeroaxis
set format x ""
unset xlabel
set grid xtics
set yrange [-1.5:1.5]
set multiplot

set size winratio,winheight

#set output scheme."-SliderCrank-v-".timestep.".pdf"
set origin 0.0,winheight*3.0+heightoff
set ylabel "corner 1"
plot resultfile  using (crankrevo($5,$7,h*$1)):((l1*sin(crankangle($5,$7,h*$1))+l2*sin(connectingangle($19,$21,h*$1))   )/c )  notitle w lines
set ylabel "corner 2"
set origin 0.0,winheight*2.0+heightoff
set ylabel "corner 3"
plot resultfile  using (crankrevo($5,$7,h*$1)):((l1*sin(crankangle($5,$7,h*$1))+l2*sin(connectingangle($19,$21,h*$1))   )/c)  notitle w lines
set ylabel "corner 4"
set origin 0.0,winheight*1.0+heightoff
plot resultfile  using (crankrevo($5,$7,h*$1)):((l1*sin(crankangle($5,$7,h*$1))+l2*sin(connectingangle($19,$21,h*$1))   )/-c) notitle w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin 0.0,winheight*0.0+heightoff
plot resultfile  using (crankrevo($5,$7,h*$1)):((l1*sin(crankangle($5,$7,h*$1))+l2*sin(connectingangle($19,$21,h*$1))   )/-c)  notitle w lines

set auto
unset ytics
unset xlabel
unset ylabel
set format

set xrange[0:2.0]
set xzeroaxis
set format x ""
set format y "%3.0e"
set multiplot
set grid xtics

set ytics axis 

set origin 0.0,winheight*3.0+heightoff
plot resultfile  using  (crankrevo($5,$7,h*$1)):($82)  t 'g1' w lines
set origin 0.0,winheight*2.0+heightoff
plot resultfile  using  (crankrevo($5,$7,h*$1)):($85)  t 'g2' w lines
set origin 0.0,winheight*1.0+heightoff
plot resultfile  using  (crankrevo($5,$7,h*$1)):($88)  t 'dot g1' w lines
set bmargin 0
set format
set xtics axis
set xlabel "Crank revolutions"
set origin 0.0,winheight*0.0+heightoff
plot resultfile  using  (crankrevo($5,$7,h*$1)):($91)  t 'dot g2' w lines


# unset multiplot
# set xzeroaxis
# set format x ""
#  unset xlabel
# set grid xtics

# set multiplot
# #set output scheme."-SliderCrank-v-".timestep.".pdf"
# set origin 0.0,winheight*3.0+heightoff
# plot resultfile  using 1:14  t 'g1' w lines
# set origin 0.0,winheight*2.0+heightoff
# plot resultfile  using 1:15  t 'g2' w lines
# set origin 0.0,winheight*1.0+heightoff
# plot resultfile  using 1:16  t 'g3' w lines
# set bmargin 0
# set format
# set xtics axis
# set xlabel "Time(s)"
# set origin 0.0,winheight*0.0+heightoff
# plot resultfile  using 1:17  t 'g4' w lines

# unset multiplot
# set xzeroaxis
# set format x ""
# unset xlabel
# set grid xtics

# set multiplot
# #set output scheme."-SliderCrank-v-".timestep.".pdf"
# set origin 0.0,winheight*3.0+heightoff
# plot resultfile  using 2:18  t 'dot g1' w lines
# set origin 0.0,winheight*2.0+heightoff
# plot resultfile  using 2:19  t 'dot g2' w lines
# set origin 0.0,winheight*1.0+heightoff
# plot resultfile  using 2:20  t 'dot g3' w lines
# set bmargin 0
# set format
# set xtics axis
# set xlabel "Crank revolutions"
# set origin 0.0,winheight*0.0+heightoff
# plot resultfile  using 2:21  t 'dot g4' w lines



# unset multiplot
# set xzeroaxis
# set format x ""
# unset xlabel
# set grid xtics

# set multiplot
# #set output scheme."-SliderCrank-v-".timestep.".pdf"
# set origin 0.0,winheight*3.0+heightoff
# plot resultfile  using 2:22  t 'lambda1' w lines
# set origin 0.0,winheight*2.0+heightoff
# plot resultfile  using 2:23  t 'lambda2' w lines
# set origin 0.0,winheight*1.0+heightoff
# plot resultfile  using 2:24  t 'lambda3' w lines
# set bmargin 0
# set format
# set xtics axis
# set xlabel "Crank revolutions"
# set origin 0.0,winheight*0.0+heightoff
# plot resultfile  using 2:25  t 'lambda4' w lines

# unset multiplot
# set xzeroaxis
# set format x ""
# unset xlabel
# set grid xtics

# set multiplot
# #set output scheme."-SliderCrank-v-".timestep.".pdf"
# set origin 0.0,winheight*3.0+heightoff
# plot resultfile  using 2:26  t 'lambda(0) 1' w lines
# set origin 0.0,winheight*2.0+heightoff
# plot resultfile  using 2:27  t 'lambda(0) 2' w lines
# set origin 0.0,winheight*1.0+heightoff
# plot resultfile  using 2:28  t 'lambda(0) 3' w lines
# set bmargin 0
# set format
# set xtics axis
# set xlabel "Crank revolutions"
# set origin 0.0,winheight*0.0+heightoff
# plot resultfile  using 2:29  t 'lambda(0) 4' w lines



# unset multiplot
# set xzeroaxis
# set format x ""
# unset xlabel
# set grid xtics

# set multiplot
# #set output scheme."-SliderCrank-v-".timestep.".pdf"
# set origin 0.0,winheight*3.0+heightoff
# plot resultfile  using 2:14  t 'g1' w lines
# set origin 0.0,winheight*2.0+heightoff
# plot resultfile  using 2:18  t 'dot g1' w lines
# set origin 0.0,winheight*1.0+heightoff
# plot resultfile  using 2:22  t 'lambda1' w lines
# set bmargin 0
# set format
# set xtics axis
# set xlabel "Crank revolutions"
# set origin 0.0,winheight*0.0+heightoff
# plot resultfile  using 2:26  t 'lambda(0) 1' w lines




