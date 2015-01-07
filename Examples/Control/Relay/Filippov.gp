#set  term X11
#set term post
#set output "Filippov.ps"
#set xrange [0:1.2]
set term tikz standalone monochrome  size 5in,5in font '\small\sf'
set output "Filippov_lambda.tex"
set xlabel "$t$"

plot\
"Filippov.dat" u 1:4 t "Siconos Platform -- INRIA $\\lambda_1$" w l,\
     "Filippov.dat" u 1:5 t "Siconos Platform -- INRIA $\\lambda_2$" w l

set output "Filippov.tex"
set xlabel "$x_1$"
set ylabel "$x_2$"
plot\
"Filippov.dat" u 2:3 t "Siconos Platform -- INRIA $x_1$ vs. $x_2$" w l


set view 82, 45, 2.0, 1.1


set zlabel "$t$"

set ticslevel 0
#set zzeroaxis  -1
set zrange [0:2.0]
# set xrange [-5.5:5.5]
# set yrange [-5.5:5.5]
set grid
#set pm3d

set output "Filippov_3D.tex"

splot\
"Filippov.dat" u 2:3:1 t "Siconos Platform -- INRIA " w l

