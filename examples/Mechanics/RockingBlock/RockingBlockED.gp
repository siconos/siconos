set autoscale
extension=".tex"

set xlabel 'time (s)' 
set grid 

resultfile = 'result.dat'
resultfile_ref = 'resultED_NewMarkAlpha.ref'


set multiplot layout 2,2
set ylabel '$y$' 
plot resultfile   using 1:3  t ' ' w lines lc rgb "red",\
     resultfile_ref   using 1:3  t ' ' w lines lc rgb "red",\

set ylabel '$\theta$'
plot resultfile    using 1:4  t ' ' w lines lc rgb "red"

set ylabel '$\dot y$'
plot resultfile   using 1:6  t ' ' w lines lc rgb "red"

set ylabel '$\dot \theta$'
plot resultfile  using 1:7   t ' ' w lines lc rgb "red"
unset multiplot

set term X11 1
set multiplot layout 2,2
set ylabel '$g_1$'
plot resultfile  using 1:8  t ' ' w lines lc rgb "red"

set ylabel '$g_2$'
plot resultfile  using 1:9  t ' ' w lines lc rgb "red"

set ylabel '$\dot g_1$'
plot resultfile  using 1:10  t ' ' w lines lc rgb "red"

set ylabel '$\dot g_2$'
plot resultfile  using 1:11  t ' ' w lines lc rgb "red"

#set ylabel '$\lambda_1$'
#plot resultfile  using 1:12  t ' ' w lines lc rgb "red"
unset multiplot

