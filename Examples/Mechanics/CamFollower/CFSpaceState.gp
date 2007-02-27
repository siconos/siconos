set  term X11
!tail -n 100000 result.dat > result-gp.dat
plot\
"result-gp.dat" u 2:3 t "Position Vs. Velocity" w l

