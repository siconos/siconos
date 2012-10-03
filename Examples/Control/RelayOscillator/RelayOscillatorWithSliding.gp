set term X11
set term pdf
extension= ".pdf"
#set term tikz standalone color solid size 5in,3in
#set term tikz standalone monochrome  size 5in,5in font '\small\sf' 
#extension=".tex"

resultfile="RelayOscillatorWithSliding.dat"

set autoscale
set ticslevel 0
set output "RelayOscillatorWithSliding-phase".extension
set view 60,240
set grid
set zlabel "x1"
set ylabel "x2"
set xlabel "x3      "

splot resultfile u ($4):3:2 w l notitle


basheight = 0.15
heightoff = 0.11
winratio = 1.0
winheight = basheight*winratio
hoffset = 0.05
voffset = 0.01
set lmargin 8 
set bmargin 0
set tmargin 0

set size 1.0 , 1.0

set xzeroaxis
set format x ""
unset xlabel
#set grid xtics
#set term tikz standalone monochrome  size 5in,3in font '\small\sf' 
set output "RelayOscillatorWithSliding".extension
set multiplot


set size winratio-0.1,winheight

#set ytics -1.0,0.5,1.
#set yrange [0:2.5]
#set xrange [0:100.0]

set origin hoffset,(winheight+voffset)*4.0+heightoff
set ylabel "x1"
plot resultfile  using 1:2  notitle w lines

set origin hoffset,(winheight+voffset)*3.0+heightoff
set ylabel "x2"
plot resultfile  using 1:3  notitle w lines

set origin hoffset,(winheight+voffset)*2.0+heightoff
set ylabel "x3"
plot resultfile  using 1:4  notitle w lines

set origin hoffset,(winheight+voffset)*1.0+heightoff
set ylabel "u"
plot resultfile  using 1:5 notitle w lines

set bmargin 0
set format
set xtics border
set xlabel "time [s]"
set ylabel "y"
set origin hoffset,winheight*0.0+heightoff

set samples 5000
plot resultfile  using 1:6  notitle w lines
unset multiplot