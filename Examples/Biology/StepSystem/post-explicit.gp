set autoscale
set term X11
#set term pdf
#extension= ".pdf"
#set term tikz standalone color solid size 5in,3in
set term tikz standalone monochrome  size 5in,5in font '\large\sf'
extension=".tex"

xmin=7.9
xmax=8.5
ymin=6.4
ymax=8.2

set xlabel "x1"
set ylabel "x2"
set xrange [xmin:xmax]
set yrange [ymin:ymax]
set output "phase-explicit-zoom1".extension
set label "\$\\theta_2^2\$" at first xmax, first 8
set label "\$\\theta_1^2\$" at first 8, first ymax+0.1

set arrow from 8,ymin to 8,ymax nohead ls 3
set arrow from xmin,8 to xmax,8 nohead ls 3

plot "result.dat" u 2:3 w l ls 1 notitle

unset label
unset arrow
xmin=7.985
xmax=8.028
ymin=7.95
ymax=8.01

set xrange [xmin:xmax]
set yrange [ymin:ymax]


set label "\$\\theta_2^2\$" at first xmax, first 8
set label "\$\\theta_1^2\$" at first 8, first ymax+0.002

set arrow from 8,ymin to 8,ymax nohead ls 3
set arrow from xmin,8 to xmax,8 nohead ls 3

set output "phase-explicit-zoom2".extension

plot "result.dat" u 2:3 w l ls 1 notitle
unset arrow
unset label

set output "phase-explicit".extension
xmin=5
xmax=10.5
ymin=5
ymax=10.5

set xrange [xmin:xmax]
set yrange [ymin:ymax]



coeff=0.05

set arrow from 8,ymin to 8,ymax nohead ls 3
set arrow from xmin,8 to xmax,8 nohead ls 3

vec1(x,y) =  ( x/(sqrt(x*x+y*y))*coeff     ) 
vec2(x,y) =  ( y/(sqrt(x*x+y*y))*coeff     ) 


set label "\$\\theta_2^2\$" at first 10.5, first 8
set label "\$\\theta_1^2\$" at first 7.7, first 10.8

plot "result.dat" u 2:3 w l ls 1 notitle

