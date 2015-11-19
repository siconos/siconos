set  term X11
!tail -n 20000 PRC.dat > PRC-gp.dat

set multiplot;                          # get into multiplot mode
set size 0.5,1;  
set origin 0.0,0.0; 
set ylabel "Vr"
set xlabel "ir"
plot "PRC-gp.dat" u 2:3 t "" w l
set origin 0.5,0.0;  
set ylabel "v0"
set xlabel "iL"
plot "PRC-gp.dat" u 4:5 t "" w l
unset multiplot                         # exit multiplot mode


