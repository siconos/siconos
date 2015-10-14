set  term X11
!tail -n 100000 result.dat > result-gp.dat
plot\
"result-gp.dat" u 1:6 t "Cam Position" w l,\
"result-gp.dat" u 1:2 t "Follower Position" w l
