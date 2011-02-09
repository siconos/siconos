set  term X11
set term post
set output "result.ps"
!tail -n 100000 resultTwoConED.dat > resultTwoConED-gp.dat
plot\
"resultTwoConED-gp.dat" u 1:2 t "Ball position" w l,\
"resultTwoConED-gp.dat" u 1:3 t "Ball Velocity" w l,\
"resultTwoConED-gp.dat" u 1:8 t "Ball acceleration" w l,\
"resultTwoConED-gp.dat" u 1:9 t "Reaction force" w l,\
"resultTwoConED-gp.dat" u 1:4 t "Reaction impulse" w l



