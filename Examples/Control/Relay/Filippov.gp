set  term X11
set term post
set output "Filippov.ps"
#set xrange [0:1.2]
plot\
"Filippov.dat" u 2:3 t "Siconos Platform -- INRIA Process x(1) vs. x(2)" w lp

plot\
"Filippov.dat" u 1:4 t "Siconos Platform -- INRIA Process lambda(1)" w lp,\
"Filippov.dat" u 1:5 t "Siconos Platform -- INRIA Process lambda(2)" w lp


