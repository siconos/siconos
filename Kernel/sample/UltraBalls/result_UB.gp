set  term X11
!tail -n 2000 resultUB.dat > resultUB-gp.dat
plot\
 "resultUB-gp.dat" u 1:2 t "Ball1 position" w l,\
 "resultUB-gp.dat" u 1:4 t "Ball2 position" w l,\
 "resultUB-gp.dat" u 1:6 t "Ball3 position" w l,\
 "resultUB-gp.dat" u 1:8 t "Ground position" w l,\
 "resultUB-gp.dat" u 1:10 t "Ceiling position" w l

# "resultUB-gp.dat" u 1:3 t "Ball1 Velocity" w l,\
# "resultUB-gp.dat" u 1:5 t "Ball2 Velocity" w l,\
# "resultUB-gp.dat" u 1:7 t "Ball3 Velocity" w l,\
# "resultUB-gp.dat" u 1:9 t "Ground Velocity" w l,\
# "resultUB-gp.dat" u 1:11 t "Ceiling Velocity" w l
