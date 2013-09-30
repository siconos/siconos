#set  term X11
plot\
"result.dat" u 1:2 t "Ball position -- q " w l,\
"result.dat" u 1:3 t "Ball Velocity -- v" w l,\
"result.dat" u 1:4 t "Impulse -- p(1)" w l,\
"result.dat" u 1:5 t "Reaction force -- p(2)" w l,\
"result.dat" u 1:6 t "Impulse multiplier -- lambda(1)" w l,\
"result.dat" u 1:7 t "Reaction multiplier -- lambda(2)" w l


