#set term X11
#tail -n 100000 result.dat > result-gp.dat
#tail -n 100000 result_tdg.dat > result_tdg-gp.dat

resultfile=  "result_tdg.dat"
resultfile1=  "result.dat"
resultfile_ref = "BouncingBallTS-D1MinusLinearOSI.ref"


plot\
 resultfile u 1:2 t "Ball position TDG" w l,\
 resultfile u 1:3 t "Ball velocity TDG" w l,\
 resultfile u 1:4 t "Reaction force TDG p(1)" w l,\
 resultfile u 1:6 t "Reaction force TDG p(2)" w l ,\
 resultfile_ref every ::1 u 1:2 t "Ball position TDG ref" w l,\
 resultfile_ref every ::1 u 1:3 t "Ball velocity TDG ref" w l,\
 resultfile_ref every ::1 u 1:4 t "Reaction force TDG p(1) ref" w l,\
 resultfile_ref every ::1 u 1:6 t "Reaction force TDG p(2) ref" w l



 # ,\
 # resultfile1 u 1:2 t "Ball position Moreau" w l,\
 # resultfile1 u 1:3 t "Ball velocity Moreau" w l,\
 # resultfile1 u 1:4 t "Reaction force Moreau" w l
 
