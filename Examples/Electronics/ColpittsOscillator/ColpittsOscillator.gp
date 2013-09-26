resultfile = 'Colpitts.dat'

resultfile = 'Colpitts.ref'

plot resultfile every ::2 u 1:2 w l t "x1", resultfile  every ::2 u 1:3 w l t "x2", resultfile   every ::2 u 1:4 w l t "x3"

set xrange [80:100]

plot resultfile  every ::2 u 1:2 w l t "x1", resultfile every ::2 u 1:3 w l t "x2", resultfile every ::2 u 1:4 w l t "x3"


plot resultfile  every ::2 u 1:5 w l t "y1", resultfile every ::2 u 1:6 w l t "y2"


plot resultfile  every ::2  u 1:7 w l t "lambda1", resultfile  every ::2 u 1:8 w l t "lambda2"
