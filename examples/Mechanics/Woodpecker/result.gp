!set term post eps color solid
!set output "woodpecker.eps"
set  term X11

set term tikz standalone monochrome  size 5in,5in font '\small\sf'
set output "woodpecker.tex"



!tail -n 100000 result.dat > result-gp.dat




set xzeroaxis
set grid xtics

set size 1.0,1.0
set origin 0.0,0.0
set multiplot

set size 0.5,0.5
set origin 0.0,0.0
plot "result-gp.dat" u 6:7 t "$\\omega_S(\\phi_S)$" w l

set size 0.5,0.5
set origin 0.5,0.0
plot "result-gp.dat" u 4:5 t "$\\omega_M(\\phi_M)$" w l

set size 0.5,0.25
set origin 0.0,0.5
plot "result-gp.dat" u 1:7 t "$\\omega_S(t)$" w l

set size 0.5,0.25
set origin 0.5,0.5
plot "result-gp.dat" u 1:5 t "$\\omega_M(t)$" w l

set size 0.5,0.25
set origin 0.0,0.75
plot "result-gp.dat" u 1:6 t "$\\phi_S(t)$" w l

set size 0.5,0.25
set origin 0.5,0.75
plot "result-gp.dat" u 1:4 t "$\\phi_M(t)$" w l

set nomultiplot
set output
