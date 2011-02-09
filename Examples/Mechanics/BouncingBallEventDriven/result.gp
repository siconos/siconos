set  term X11
plot\
"result.dat" u 1:2 t "Ball position" w l,\
"result.dat" u 1:3 t "Ball Velocity" w l,\
"result.dat" u 1:4 t "Impulse" w l,\
"result.dat" u 1:5 t "Reaction force" w l


