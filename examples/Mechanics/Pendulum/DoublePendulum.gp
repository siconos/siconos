set  term X11
set multiplot


basheight = 0.45
heightoff = -0.85
winratio = 1.0
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

set origin 0.0,winheight*3.0+heightoff
 

#plot\
#"DoublePendulumResult.dat" u 1:2 t "Mass 1 position" w l,\
#"DoublePendulumResult.dat" u 1:3 t "Mass 1 Velocity" w l,\
#"DoublePendulumResult.dat" u 1:4 t "Mass 2 position" w l,\
#"DoublePendulumResult.dat" u 1:5 t "Mass 2 Velocity" w l,\
#"DoublePendulumResult.dat" u 1:6 t "x1" w l,\
#"DoublePendulumResult.dat" u 1:8 t "x2" w l

plot\
"DoublePendulumResult.dat" u 1:10 t "v1" w l,\
"DoublePendulumResult.dat" u 1:11 t "v2" w l

set size 0.4 , 0.4
set origin 0.0,winheight*2.0+heightoff-0.01

set xrange [-1:3]
set yrange [-2:2]

plot\
"DoublePendulumResult.dat" u 6:7 t "Mass 1 trajectory" w l,\
"DoublePendulumResult.dat" u 8:9 t "Mass 2 trajectory" w l


