set  term X11

plot\
 "ObserverLCS.dat" u 1:2 t "x(1)" w l,\
 "ObserverLCS.dat" u 1:3 t "x(2)" w l,\
 "ObserverLCS.dat" u 1:4 t "hatx(1)" w l,\
 "ObserverLCS.dat" u 1:5 t "hatx(2)" w l,\
 "ObserverLCS.dat" u 1:6 t "lambda" w l,\
 "ObserverLCS.dat" u 1:7 t "hatlambda" w l
 

