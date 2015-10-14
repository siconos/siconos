set autoscale

# choose terminal type
aqua=1
X11=0

set xlabel 'time (s)' 
set grid 

data = 'result.dat'
dataref = 'result_AlphaScheme.ref'
color = "blue"
colorref = "red"
method = "result"


#set terminal pdf
#set output "figure1.pdf"

if (X11==1) set term X11 0 
if (aqua==1) set term aqua 0
set multiplot layout 2,2
set ylabel '$x$' 

plot \
     data using 1:2 t method w lines lc rgb color,\
     dataref every ::2 using 1:2 t 'ref' w lines lc rgb colorref
     

set ylabel '$y$'
plot data using 1:3 t method w lines lc rgb color,\
     dataref every ::2 using 1:3 t 'ref' w lines lc rgb colorref


set ylabel '$\dot x$'
plot data using 1:4 t method w lines lc rgb color,\
     dataref every ::2 using 1:4 t 'ref' w lines lc rgb colorref


set ylabel '$\dot y$'
plot data using 1:5 t method w lines lc rgb color,\
     dataref every ::2 using 1:5 t 'ref' w lines lc rgb colorref
unset multiplot

#set terminal pdf
#set output "figure2.pdf"

if (X11==1)set term X11 1 
if (aqua==1) set term aqua 1

set multiplot layout 2,2
set ylabel '$g$'
plot data using 1:8 t method w lines lc rgb color


set ylabel '$\dot g$'
plot data using 1:9 t method w lines lc rgb color


set ylabel '$\lambda$'
plot data using 1:10 t method w lines lc rgb color
unset multiplot
