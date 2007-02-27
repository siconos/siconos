set  term X11
!tail -n 2000 result.dat > result-gp.dat
plot\
 "result-gp.dat" u 1:2 t "Bar position x" w l,\
 "result-gp.dat" u 1:3 t "Bar position y" w l, \
 "result-gp.dat" u 1:4 t "Bar rotation theta" w l

# "result-gp.dat" u 1:4 t "Ground position" w l,\
# "result-gp.dat" u 1:5 t "Ground Velocity" w l, \


