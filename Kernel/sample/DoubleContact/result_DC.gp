set  term X11
!tail -n 2000 resultDC.dat > resultDC-gp.dat
plot\
 "resultDC-gp.dat" u 1:2 t "Ball1 position" w l,\
 "resultDC-gp.dat" u 1:3 t "Ball1 Velocity" w l,\
 "resultDC-gp.dat" u 1:4 t "Ball2 position" w l,\
 "resultDC-gp.dat" u 1:5 t "Ball2 Velocity" w l,\
 "resultDC-gp.dat" u 1:6 t "Ball3 position" w l,\
 "resultDC-gp.dat" u 1:7 t "Ball3 Velocity" w l
 
# "resultDC-gp.dat" u 1:6 t "Reaction force" w l

# "resultDC-gp.dat" u 1:4 t "Ground position" w l,\
# "resultDC-gp.dat" u 1:5 t "Ground Velocity" w l, \


