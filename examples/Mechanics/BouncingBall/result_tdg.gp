#set term X11
#tail -n 100000 result.dat > result-gp.dat
#tail -n 100000 result_tdg.dat > result_tdg-gp.dat
plot\
"result_tdg-gp.dat" u 1:2 t "Ball position TDG" w l,\
"result_tdg-gp.dat" u 1:3 t "Ball velocity TDG" w l,\
"result_tdg-gp.dat" u 1:4 t "Reaction force TDG p(1)" w l,\
"result_tdg-gp.dat" u 1:6 t "Reaction force TDG p(2)" w l,\
"result-gp.dat" u 1:2 t "Ball position Moreau" w l,\
"result-gp.dat" u 1:3 t "Ball velocity Moreau" w l,\
"result-gp.dat" u 1:4 t "Reaction force Moreau" w l
