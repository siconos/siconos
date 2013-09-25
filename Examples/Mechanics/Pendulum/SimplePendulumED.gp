
set autoscale

set xlabel 'time (s)' 
set grid 

data = 'result.dat'
#data = 'result_AlphaScheme.ref'
color = "blue"
method = " "
#set terminal pdf
#set output "figure1.pdf"
#set term X11 0 
set multiplot layout 2,2
set ylabel '$x$' 

plot data using 1:2 t method w lines lc rgb color

set ylabel '$y$'
plot data using 1:3 t method w lines lc rgb color


set ylabel '$\dot x$'
plot data using 1:4 t method w lines lc rgb color


set ylabel '$\dot y$'
plot data using 1:5 t method w lines lc rgb color
unset multiplot

#set terminal pdf
#set output "figure2.pdf"
#set term X11 1
set multiplot layout 2,2
set ylabel '$g$'
plot data using 1:8 t method w lines lc rgb color


set ylabel '$\dot g$'
plot data using 1:9 t method w lines lc rgb color


set ylabel '$\lambda$'
plot data using 1:10 t method w lines lc rgb color
unset multiplot
