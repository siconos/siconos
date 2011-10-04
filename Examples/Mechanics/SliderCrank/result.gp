set term X11
!tail -n 100000 result.dat > result-gp.dat
#plot "result-gp.dat" u 1:2 t "Crank revolutions" w l
plot "result-gp.dat" u 2:3 t "Rod angle" w l
#plot "result-gp.dat" u 2:4 t "Slider angle" w l
#plot "result-gp.dat" u 2:5 t "Crank angular velocity" w l
#plot "result-gp.dat" u 2:6 t "Rod angular velocity" w l
#plot "result-gp.dat" u 2:7 t "Slider angular velocity" w l
#plot "result-gp.dat" u 2:8 t "Corner 1 y-position (normalized)" w l
#plot "result-gp.dat" u 2:9 t "Corner 2 y-position (normalized)" w l
#plot "result-gp.dat" u 2:10 t "Corner 3 y-position (normalized)" w l
#plot "result-gp.dat" u 2:11 t "Corner 4 y-position (normalized)" w l
#plot "result-gp.dat" u 2:12 t "Slider x-position (normalized)" w l
#plot "result-gp.dat" u 2:13 t "Slider y-position (normalized)" w l
#plot "result-gp.dat" u 3:6 t "Phase diagram rod" w l
#plot "result-gp.dat" u 12:13 t "Diagram slider position" w l
