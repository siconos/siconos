set  term X11
!tail -n 100000 result.dat > result-gp.dat
plot\
"result-gp.dat" u 1:2 t "Box position" w l,\
"result-gp.dat" u 1:3 t "Box Velocity" w l,\
"result-gp.dat" u 1:4 t "Reaction Force 1" w l,\
"result-gp.dat" u 1:5 t "Reaction Force 2" w l,\
"result-gp.dat" u 1:6 t "Reaction Force 3" w l,\
"result-gp.dat" u 1:7 t "Reaction Force 4" w l
