set term X11
!tail -n 100000 result.dat > result-gp.dat
#plot "result-gp.dat" u 1:8 t "Corner 1" w l
#plot "result-gp.dat" u 1:9 t "Corner 2" w l
#plot "result-gp.dat" u 1:10 t "Corner 3" w l
#plot "result-gp.dat" u 1:11 t "Corner 4" w l
plot "result-gp.dat" u 1:5 t "Crank speed" w l
#plot "result-gp.dat" u 1:6 t "Rod speed" w l
#plot "result-gp.dat" u 3:6 t "Phase diagram rod" w l
#plot "result-gp.dat" u 12:13 t "Phase diagram slider position" w l
