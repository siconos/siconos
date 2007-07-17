set  term X11
!tail -n 100000 result.dat > result-gp.dat
plot\
"result-gp.dat" u 4:2 t "Ball position" w l
