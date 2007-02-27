set  term X11
!tail -n 20000 result.dat > result-gp.dat
set lmargin 2
set size 1,0.25
set multiplot
set origin 0.0,0.75
plot\
 "result-gp.dat" u 1:2 t "Ball1 position z" w l,\
 "result-gp.dat" u 1:3 t "Ball2 position z w l
set origin 0.0,0.50
plot  "result-gp.dat" u 1:2 t "Ball1 Velocity" w l,\
 "result-gp.dat" u 1:3 t "Ball2 Velocity" w l
set nomultiplot


