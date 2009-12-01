set  term X11
set term post
set output "OptimalControl.ps"
#set grid
#set xrange [0:1.2]
plot\
"OptimalControl.dat" u 1:2 t "Siconos Platform -- INRIA                                             Process x(1)" w l,\
"OptimalControl.dat" u 1:3 t "Process x(2)" w l,\
"OptimalControl.dat" u 1:4 t "Process p(1)" w l,\
"OptimalControl.dat" u 1:5 t "Process p(2)" w l

plot\
"OptimalControl.dat" u 1:6 t "lambda(1)" w l,\
"OptimalControl.dat" u 1:7 t "lambda(2)" w l,\
"OptimalControl.dat" u 1:8 t "y(1)" w l,\
"OptimalControl.dat" u 1:9 t "y(2)" w l


plot\
"OptimalControl.dat" u 1:(0.5*($6-$9)) t "lambda" w lp,\
"OptimalControl.dat" u 1:($8-$7) t "y" w lp
