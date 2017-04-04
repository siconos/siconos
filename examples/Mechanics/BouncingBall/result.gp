#set  term X11
#!tail -n 100000 result.dat > result-gp.dat
resultfile = "result.dat"
#resultfile = "BouncingBallSchatzmanPaoliOSI.ref"
#resultfile = "result_tdg.dat"
plot \
resultfile every ::2 u 1:2 t "Ball position" w l,\
resultfile every ::2 u 1:3 t "Ball Velocity" w l,\
resultfile every ::2 u 1:4 t "Reaction force" w l

#,\
#"result-gp.dat" u 1:5 t "Multiplier" w l
     


