set  term X11
!tail -n 100000 result.dat > result-gp.dat
set lmargin 10 
set size 1,0.25
set multiplot
set origin 0.0,0.75
plot\
"result-gp.dat" u 1:2 t "Ball position" w l,\
"result-gp.dat" u 1:3 t "Ball Velocity" w l,\
"result-gp.dat" u 1:4 t "Reaction force" w l
set origin 0.0,0.50
plot\
"result-gp.dat" u 1:5 t "Ball position" w l,\
"result-gp.dat" u 1:6 t "Ball Velocity" w l,\
"result-gp.dat" u 1:7 t "Reaction force" w l
set nomultiplot
