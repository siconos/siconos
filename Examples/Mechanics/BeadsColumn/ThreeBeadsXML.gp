set  term X11
!tail -n 20000 result.dat > result-gp.dat
set lmargin 10 
set size 1,0.25
set multiplot
set origin 0.0,0.75
plot\
 "result-gp.dat" u 1:2 t "Ball1 position" w l,\
 "result-gp.dat" u 1:4 t "Ball2 position" w l,\
 "result-gp.dat" u 1:6 t "Ball3 position" w l
set origin 0.0,0.50
 plot "result-gp.dat" u 1:3 t "Ball1 Velocity" w l,\
 "result-gp.dat" u 1:5 t "Ball2 Velocity" w l,\
 "result-gp.dat" u 1:7 t "Ball3 Velocity" w l
set origin 0.0,0.25
plot  "result-gp.dat" u 1:2 t "Ball1 position" w l,\
 "result-gp.dat" u 1:3 t "Ball1 Velocity" w l
set origin 0.0,0.0
plot  "result-gp.dat" u 1:4 t "Ball2 position" w l,\
 "result-gp.dat" u 1:5 t "Ball2 Velocity" w l
set nomultiplot
