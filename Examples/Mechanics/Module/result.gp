set  term X11

plot\
 "result.dat" u 1:2 t " Bar reaction force" w l, \
 "result.dat" u 1:3 t " Bar reaction force" w l, \
 "result.dat" u 1:4 t " Bar reaction force" w l, \
 "result.dat" u 1:5 t " Bar reaction force" w l, \
 "result.dat" u 1:6 t " Bar reaction force" w l, \
 "result.dat" u 1:7 t " Bar reaction force" w l, \
 "result.dat" u 1:8 t " Bar reaction force" w l, \
 "result.dat" u 1:9 t " Bar reaction force" w l


 # "result-gp.dat" u 1:5 t "Bar velocity x" w l,\
 #"result-gp.dat" u 1:6 t "velocity y" w l, \
 # "result-gp.dat" u 1:7 t "velocity theta" w l,\
#!tail -n 2000 result.dat > result-gp.dat

pause -1
set term postscript eps color
set output "result.eps"
replot

