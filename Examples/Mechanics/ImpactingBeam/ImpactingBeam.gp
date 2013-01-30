set  term X11

set term pdf

set output "ImpactingBeam-Energy.pdf"
set auto
plot\
"ImpactingBeam.dat" u 1:6 t "Potential energy" w l,\
"ImpactingBeam.dat" u 1:7 t "Kinetic energy" w l,\
"ImpactingBeam.dat" u 1:($7+$6) t "Total energy" w l

set output "ImpactingBeam-qvr.pdf"


basheight = 0.3
heightoff = 0.05
winratio = 1.0
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
set origin 0.01,winheight*2.0+heightoff
plot\
"ImpactingBeam.dat" u 1:2 t "Beam end position" w l


set ylabel "m/s" 
set origin 0.01,winheight+heightoff
plot\
"ImpactingBeam.dat" u 1:3 t "Beam end Velocity" w l

set bmargin 0
set format
set xtics axis
set xlabel "time in s"

set ylabel "Ns" 
set origin 0.01,0.0+heightoff

plot\
"ImpactingBeam.dat" u 1:4 t "Reaction force" w l

set nomultiplot
set origin


