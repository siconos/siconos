#set  term X11
set multiplot

resultfile = 'SimplePendulumResult.dat'
resultfile_ref = 'result_AlphaScheme.ref'

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
#"SimplePendulumResult.dat" u 1:2 t "Mass 1 position" w l,\
#"SimplePendulumResult.dat" u 1:3 t "Mass 1 Velocity" w l,\
#"SimplePendulumResult.dat" u 1:4 t "x1" w l

plot \
resultfile u 1:6 t "v1" w l,\
resultfile_ref every ::1 u 1:6 t "v1_ref" w l




set origin 0.0,winheight*2.0+heightoff-0.01

set xrange [-0.0:2]
set yrange [-2:2]

plot resultfile u 4:5 t "Mass 1 trajectory" w l,\
     resultfile_ref every ::1  u 4:5 t "Mass 1 trajectory (ref)" w l,\

