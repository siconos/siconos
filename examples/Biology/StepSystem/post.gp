set autoscale
set term X11
#set term pdf
#extension= ".pdf"
#set term tikz standalone color solid size 5in,3in
set term tikz standalone monochrome  size 5in,5in font '\small\sf' 
extension=".tex"


resultfile = "simu.1.6.log"
#resultfile = "simu.1.6.ref"


set autoscale
set xrange [0:2.0]
set output "traj".extension
plot resultfile u 1:2 w l ,\
     resultfile u 1:3 w l ,\
     resultfile u 1:8 w l ,\
     resultfile u 1:9 w l 
     

set output "phase".extension

set xrange [0:12.5]
set yrange [0:12.5]
coeff=0.05

set arrow from 4,0 to 4,12.5 nohead ls 3
set arrow from 8,0 to 8,12.5 nohead ls 3
set arrow from 0,4 to 12.5,4 nohead ls 3
set arrow from 0,8 to 12.5,8 nohead ls 3

vec1(x,y) =  ( x/(sqrt(x*x+y*y))*coeff     ) 
vec2(x,y) =  ( y/(sqrt(x*x+y*y))*coeff     ) 

set xlabel "\$x_1\$"
set ylabel "\$x_2\$"

set label "\$\\theta_2^1\$" at first 12.5, first 4
set label "\$\\theta_2^2\$" at first 12.5, first 8
set label "\$\\theta_1^1\$" at first 3.7, first 12.8
set label "\$\\theta_1^2\$" at first 7.7, first 12.8
plot "simu.1.6.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.1.6.log" u 2:3 w l ls 1 notitle,\
     "simu.1.7.5.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.1.7.5.log" u 2:3 w l ls 1 notitle,\
     "simu.1.9.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.1.9.log" u 2:3 w l ls 1 notitle,\
     "simu.1.11.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.1.11.log" u 2:3 w l ls 1 notitle,\
     "simu.2.4.5.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.2.4.5.log" u 2:3 w l ls 1 notitle,\
     "simu.3.1.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.3.1.log" u 2:3 w l ls 1 notitle,\
     "simu.3.11.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.3.11.log" u 2:3 w l ls 1 notitle,\
     "simu.3.5.0.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.3.5.0.log" u 2:3 w l ls 1 notitle,\
     "simu.4.4.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.4.4.log" u 2:3 w l ls 1 notitle,\
     "simu.4.1.4.1.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.4.1.4.1.log" u 2:3 w l ls 1 notitle,\
     "simu.5.0.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.5.0.log" u 2:3 w l ls 1 notitle,\
     "simu.6.5.10.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.6.5.10.log" u 2:3 w l ls 1 notitle,\
     "simu.7.12.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.7.12.log" u 2:3 w l ls 1 notitle,\
     "simu.7.1.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.7.1.log" u 2:3 w l ls 1 notitle,\
     "simu.10.11.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.10.11.log" u 2:3 w l ls 1 notitle,\
     "simu.10.1.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.10.1.log" u 2:3 w l ls 1 notitle,\
     "simu.11.10.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.11.10.log" u 2:3 w l ls 1 notitle,\
     "simu.11.5.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.11.5.log" u 2:3 w l ls 1 notitle,\
     "simu.11.3.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.11.3.log" u 2:3 w l ls 1 notitle,\
     "simu.12.9.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.12.9.log" u 2:3 w l ls 1 notitle,\
     "simu.12.8.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.12.8.log" u 2:3 w l ls 1 notitle,\
     "simu.12.6.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "simu.12.6.log" u 2:3 w l ls 1 notitle
set output "phase-reduced".extension




# plot  "simu.11.8.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
#       "simu.11.8.log" u 2:3 w l ls 1 notitle





