set  term X11
!tail -n 2000 resultB2D.dat > resultB2D-gp.dat
plot\
 "resultB2D-gp.dat" u 3:2 t "Ball position" w l,\
 "resultB2D-gp.dat" u 1:4 t "Ground position" w l
 
# "resultB2D-gp.dat" u 1:2 t "Ball position x" w l,\
# "resultB2D-gp.dat" u 1:3 t "Ball position y" w l,\
# "resultB2D-gp.dat" u 1:4 t "Ground position" w l
 
# "resultB2D-gp.dat" u 1:6 t "Reaction force" w l
# "resultB2D-gp.dat" u 1:4 t "Ground position" w l,\
# "resultB2D-gp.dat" u 1:5 t "Ground Velocity" w l, \


