#set  term X11
!tail -n 100000 result.dat > result-gp.dat
plot\
"result-gp.dat" u 7:4 t "Position 1" w l


