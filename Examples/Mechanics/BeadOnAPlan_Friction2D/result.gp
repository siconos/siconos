set  term X11
!tail -n 100000 result.dat > result-gp.dat
plot\
"result-gp.dat" u 1:2 t "position" w l,\
"result-gp.dat" u 1:3 t "Velocity" w l