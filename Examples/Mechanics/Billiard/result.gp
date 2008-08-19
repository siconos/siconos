set  term X11
!tail -n 20000 result.dat > result-gp.dat

plot\
"result-gp.dat" u 1:2 t "Ball position1" w l,\
"result-gp.dat" u 1:3 t "Ball Velocity1" w l,\
"result-gp.dat" u 1:4 t "Ball reaction1" w l,\
"result-gp.dat" u 1:5 t "Ball position2" w l,\
"result-gp.dat" u 1:6 t "Ball Velocity2" w l,\
"result-gp.dat" u 1:7 t "Ball reaction2" w l