set autoscale
set term X11
#set term pdf
#extension= ".pdf"
#set term tikz standalone color solid size 5in,3in
set term tikz standalone monochrome  size 5in,5in font '\small\sf' 
extension=".tex"






set autoscale
set xrange [0:0.25]
set output "traj-sign".extension
plot "sign.11.8.log" u 1:2 w l ,\
     "sign.11.8.log" u 1:3 w l ,\
     "sign.11.8.log" u 1:8 w l ,\
     "sign.11.8.log" u 1:9 w l 
     

set output "phase-sign".extension

set xrange [0:12.5]
set yrange [0:12.5]
coeff=0.05

set arrow from 4,0 to 4,12.5 nohead ls 3
set arrow from 8,0 to 8,12.5 nohead ls 3
set arrow from 0,4 to 12.5,4 nohead ls 3
set arrow from 0,8 to 12.5,8 nohead ls 3

vec1(x,y) =  ( x/(sqrt(x*x+y*y))*coeff     ) 
vec2(x,y) =  ( y/(sqrt(x*x+y*y))*coeff     ) 
plot "sign.1.6.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.1.6.log" u 2:3 w l ls 1 notitle,\
     "sign.1.7.5.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.1.7.5.log" u 2:3 w l ls 1 notitle,\
     "sign.1.9.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.1.9.log" u 2:3 w l ls 1 notitle,\
     "sign.1.11.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.1.11.log" u 2:3 w l ls 1 notitle,\
     "sign.2.4.5.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.2.4.5.log" u 2:3 w l ls 1 notitle,\
     "sign.3.1.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.3.1.log" u 2:3 w l ls 1 notitle,\
     "sign.3.11.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.3.11.log" u 2:3 w l ls 1 notitle,\
     "sign.3.5.0.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.3.5.0.log" u 2:3 w l ls 1 notitle,\
     "sign.4.4.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.4.4.log" u 2:3 w l ls 1 notitle,\
     "sign.5.0.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.5.0.log" u 2:3 w l ls 1 notitle,\
     "sign.6.5.10.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.6.5.10.log" u 2:3 w l ls 1 notitle,\
     "sign.7.12.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.7.12.log" u 2:3 w l ls 1 notitle,\
     "sign.7.1.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.7.1.log" u 2:3 w l ls 1 notitle,\
     "sign.10.11.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.10.11.log" u 2:3 w l ls 1 notitle,\
     "sign.10.1.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.10.1.log" u 2:3 w l ls 1 notitle,\
     "sign.11.10.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.11.10.log" u 2:3 w l ls 1 notitle,\
     "sign.11.5.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.11.5.log" u 2:3 w l ls 1 notitle,\
     "sign.11.3.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.11.3.log" u 2:3 w l ls 1 notitle,\
     "sign.12.9.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.12.9.log" u 2:3 w l ls 1 notitle,\
     "sign.12.8.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.12.8.log" u 2:3 w l ls 1 notitle,\
     "sign.12.6.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
     "sign.12.6.log" u 2:3 w l ls 1 notitle


# plot "sign.11.8.log" u 2:3:(vec1($8,$9)):(vec2($8,$9)) w vec notitle,\
#      "sign.11.8.log" u 2:3 w l notitle
    







