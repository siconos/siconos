set  term X11
!tail -n 2000 resultRB.dat > resultRB-gp.dat
plot\
 "resultRB-gp.dat" u 1:2 t "Ball1 position" w l,\
 "resultRB-gp.dat" u 1:3 t "Ball1 Velocity" w l,\
 "resultRB-gp.dat" u 1:4 t "Ball2 position" w l,\
 "resultRB-gp.dat" u 1:5 t "Ball2 Velocity" w l,\
 "resultRB-gp.dat" u 1:6 t "Ball3 position" w l,\
 "resultRB-gp.dat" u 1:7 t "Ball3 Velocity" w l
 
# "resultRB-gp.dat" u 1:6 t "Reaction force" w l

# "resultRB-gp.dat" u 1:4 t "Ground position" w l,\
# "resultRB-gp.dat" u 1:5 t "Ground Velocity" w l, \


