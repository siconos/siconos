set  term X11
!tail -n 2000 resultSF.dat > resultSF-gp.dat
plot\
 "resultSF-gp.dat" u 3:2 t "Ball position" w l,\
 "resultSF-gp.dat" u 1:4 t "Ground position" w l
