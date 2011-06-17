set  term X11
set term post
set output "SimpleExampleRelay.ps"
set xrange [0:1.2]
plot\
"SimpleExampleRelay.dat" u 1:2 t "Siconos Platform -- INRIA                                             Process x(1)" w l,\
"SimpleExampleRelay.dat" u 1:3 t "Process x(2)" w l,\
"SimpleExampleRelay.dat" u 1:6 t "lambda(1)" w l,\
"SimpleExampleRelay.dat" u 1:7 t "lambda(2)" w l


