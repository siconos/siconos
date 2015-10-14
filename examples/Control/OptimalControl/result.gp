set  term X11
set term post
set output "OptimalControl.ps"
set grid
#set xrange [0:1.2]

plot\
"OptimalControl.dat" u 1:2 t "Siconos Platform -- INRIA                                             Process x(1)" w l linecolor rgb "cyan",\
"OptimalControl.dat" u 1:3 t "Process x(2)" w l linecolor rgb "blue" ,\
"trajectoires.dat" u 5:1 t "x(1)_arc" w l linecolor rgb "red",\
"trajectoires.dat" u 5:2 t "x(2)__arc" w l linecolor rgb "green"

plot\
"OptimalControl.dat" u 1:4 t "Process p(1)" w l linecolor rgb "cyan",\
"OptimalControl.dat" u 1:5 t "Process p(2)" w l linecolor rgb "blue" ,\
"trajectoires.dat" u 5:3 t "p(1)_arc" w l linecolor rgb "red",\
"trajectoires.dat" u 5:4 t "p(2)__arc" w l linecolor rgb "green"

plot\
"OptimalControl.dat" u 1:6 t "lambda(1)" w l linecolor rgb "cyan",\
"OptimalControl.dat" u 1:7 t "lambda(2)" w l linecolor rgb "blue",\
"OptimalControl.dat" u 1:8 t "y(1)" w l linecolor rgb "red",\
"OptimalControl.dat" u 1:9 t "y(2)" w l linecolor rgb "green"


plot\
"OptimalControl.dat" u 1:(0.5*($6-$9)) t "lambda" w l linecolor rgb "green",\
"OptimalControl.dat" u 1:($8-$7) t "y_siconos" w l linecolor rgb "red",\
"trajectoires.dat"   u 5:6 t "y_traj_sing" w l linecolor rgb "blue"

plot\
"OptimalControl.dat" u 2:3 t "x(2) = f(x(1)) siconos " w l, \
"trajectoires.dat"   u 1:2 t "x(2) = f(x(1)) traj_singuls " w l

set yrange [-2:2]
plot\
"OptimalControl.dat" u 1:(0.5*(1.0+(1.0-$9))) t "u(t)" w l linecolor rgb "red"


