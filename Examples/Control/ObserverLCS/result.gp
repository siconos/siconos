set  term X11
!tail -n 100000 ObserverLCS.dat > result-gp.dat
set xrange [0:1]
plot\
"result-gp.dat" u 1:2 t "Process x(1)" w l,\
"result-gp.dat" u 1:3 t "Process x(2)" w l,\
"result-gp.dat" u 1:4 t "Observer x(1)" with l ,\
"result-gp.dat" u 1:5 t "Observer x(2)" w l,\
"result-gp.dat" u 1:6 t "lambda" w l,\
"result-gp.dat" u 1:7 t "hatLambda" w l


