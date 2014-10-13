set  term X11
set term pdf
set term pdf monochrome dashed size 5in,3.5in
extension= ".pdf"
set term tikz standalone color solid size 5in,3in
set term tikz standalone monochrome  size 5in,3in font '\small\sf'
extension=".tex"


resultfile='ImpactingBar.dat'
method="MoreauJean-theta05-"

method="D1MinusLinear-"
resultfile='ImpactingBarD1MinusLinear.dat'

model="ImpactingBar-damping-"
timestep = "1e-7"
h=1e-7
ndof = "1000-"


file = model.method.ndof.timestep


set output file."-Energy".extension
set auto
plot \
     resultfile u 1:6 t "Potential energy" w l,\
     resultfile u 1:7 t "Kinetic energy" w l,\
     resultfile u 1:($7+$6) t "Total energy" w l


set output file."-Caplambda".extension
plot \
     resultfile u 1:4 t "Reaction impulse" w l



set output file."-lambda".extension
plot \
      resultfile u 1:($12*2/h) t "Reaction force (lambda)" w l ,\
      resultfile u 1:($13*2/h) t "Reaction force (lambda)" w l



set output file."-qvr".extension
offset=0.05
xoff=0.05
set format y "%.1e"

basheight = 0.3
heightoff = 0.07
winratio = 0.95
winheight = basheight*winratio

set lmargin 8 
set bmargin 0
set tmargin 0

set size 1.0 , 1.0

set xzeroaxis
set format x ""

set mxtics 
set grid xtics
set grid mxtics

set multiplot
set size winratio,winheight

set ylabel "m" 
set origin xoff,winheight*2.0+heightoff+offset+offset
plot \
     resultfile u (1000*$1):2 t "q" w l
set ylabel "m/s"
set origin xoff,winheight+heightoff+offset
plot \
     resultfile u (1000*$1):3 t "v" w l
set format x "%.2f"
#set bmargin 0
set xtics border
set xlabel "time [ms]"
set ylabel "N"

set origin xoff,0.0+heightoff
plot \
     resultfile u (1000*$1):($4) t "r" w l

set nomultiplot
set output file."-qvrend".extension
set multiplot
set key right bottom

set xzeroaxis
set format x ""
set xlabel ""

set ylabel "m"
set origin xoff,winheight*2.0+heightoff+offset+offset
plot \
     resultfile u (1000*$1):8 t "q (end of the bar)" w l
set ylabel "m/s" 
set origin xoff,winheight+heightoff+offset
plot \
     resultfile u (1000*$1):9 t "v (end of the bar)" w l
set format x "%.2f"
#set bmargin 0
set xtics border
set xlabel "time [ms]"
set ylabel "N"

set origin xoff,0.0+heightoff
plot \
     resultfile u (1000*$1):($11) t "v (middle of the bar)" w l





