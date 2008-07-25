set  term X11
set term post
set output "SingleDSObserverLCS-jump.ps"
set xrange [-0.1:30]
#set xrange [-0.1:3]


plot\
 "SingleDSObserverLCS.dat" u 1:4 t "Siconos Platform -- INRIA                                              Process x(1)" w l,\
 "SingleDSObserverLCS.dat" u 1:5 t "Process x(2)" w l,\
 "SingleDSObserverLCS.dat" u 1:2 t "Observer x(1)" w l,\
 "SingleDSObserverLCS.dat" u 1:3 t "Observer x(2)" w l,\
 "SingleDSObserverLCS.dat" u 1:8 t "u" w l
set output "SingleDSObserverLCS-jump-1.ps"
 # "SingleDSObserverLCS.dat" u 1:6 t "Observer x(1)" w l,\
# "SingleDSObserverLCS.dat" u 1:7 t "Observer x(2)" w l,\

set xrange [-0.1:25]
plot\
 "SingleDSObserverLCS.dat" u 1:9 t "Siconos Platform -- INRIA                                  Observation Error x(1)" w l,\
  "SingleDSObserverLCS.dat" u 1:10 t "                                                         Observation Error x(2)" w l
set output "SingleDSObserverLCS-jump-2.ps"
set xrange [20:26]
set yrange [-0.3:0.3]
plot\
 "SingleDSObserverLCS.dat" u 1:4 t "Siconos Platform -- INRIA                                              Process x(1)" w l,\
 "SingleDSObserverLCS.dat" u 1:5 t "Process x(2)" w l,\
 "SingleDSObserverLCS.dat" u 1:2 t "Observer x(1)" w l,\
 "SingleDSObserverLCS.dat" u 1:3 t "Observer x(2)" w l,\
 "SingleDSObserverLCS.dat" u 1:8 t "u" w l
