set  term X11
#set term post
set output "SingleDSObserverLCS.ps"
set xrange [0:1.2]
plot\
 "SingleDSObserverLCS.dat" u 1:4 t "Process x(1)" w l,\
 "SingleDSObserverLCS.dat" u 1:5 t "Process x(2)" w l,\
 "SingleDSObserverLCS.dat" u 1:2 t "Observer x(1)" w l,\
 "SingleDSObserverLCS.dat" u 1:3 t "Observer x(2)" w l,\
 "SingleDSObserverLCS.dat" u 1:6 t "lambda" w l,\
 "SingleDSObserverLCS.dat" u 1:7 t "hatlambda" w l
 

