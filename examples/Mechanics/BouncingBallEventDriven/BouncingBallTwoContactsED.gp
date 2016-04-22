set term X11
set term pdf
set output "BouncingBallTwoContactsED.pdf"
resultfile = "BouncingBallTwoContactsED.dat"
#resultfile = "BouncingBallTwoContactsED.ref"
plot \
resultfile u 1:2 every ::1 t "Ball position" w l,\
resultfile u 1:3 every ::1 t "Ball Velocity" w l,\
resultfile u 1:8 every ::1 t "Ball acceleration" w l,\
resultfile u 1:9 every ::1 t "Reacti on force" w l,\
resultfile u 1:4 every ::1 t "Reaction impulse" w l



