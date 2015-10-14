set autoscale
extension=".tex"

set xlabel 'time (s)' 
set grid 

set term X11 0
set multiplot layout 2,2
set ylabel '$y$' 
plot 'result.dat'   using 1:3  t ' ' w lines lc rgb "red"

set ylabel '$\theta$'
plot 'result.dat'    using 1:4  t ' ' w lines lc rgb "red"

set ylabel '$\dot y$'
plot 'result.dat'   using 1:6  t ' ' w lines lc rgb "red"

set ylabel '$\dot \theta$'
plot 'result.dat'  using 1:7   t ' ' w lines lc rgb "red"
unset multiplot

set term X11 1
set multiplot layout 2,2
set ylabel '$g_1$'
plot 'result.dat'  using 1:8  t ' ' w lines lc rgb "red"

set ylabel '$g_2$'
plot 'result.dat'  using 1:9  t ' ' w lines lc rgb "red"

set ylabel '$\dot g_1$'
plot 'result.dat'  using 1:10  t ' ' w lines lc rgb "red"

set ylabel '$\dot g_2$'
plot 'result.dat'  using 1:11  t ' ' w lines lc rgb "red"

#set ylabel '$\lambda_1$'
#plot 'result.dat'  using 1:12  t ' ' w lines lc rgb "red"
unset multiplot

