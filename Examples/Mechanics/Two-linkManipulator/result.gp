set  term X11
!tail -n 100000 result.dat > result-gp.dat
plot\
"result-gp.dat" u 1:2 t "Position 1" w l,\
"result-gp.dat" u 1:3 t "Velocity 1" w l,\
"result-gp.dat" u 1:4 t "Position 2" w l,\
"result-gp.dat" u 1:5 t "Velocity 2" w l,\
"result-gp.dat" u 1:6 t "y" w l


