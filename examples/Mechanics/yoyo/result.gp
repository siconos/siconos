
set  term X11
!tail -n 100000 result.dat > result-gp.dat


set size 1,1./3.
set multiplot


set origin 0.0,2./3.
plot\
"result-gp.dat" u 1:2 t "Yoyo position teta" w l,\
"result-gp.dat" u 1:3 t "Yoyo vitesse teta" w l


set origin 0.0,1./3.
plot\
"result-gp.dat" u 1:6 t "Yoyo position y" w l,\
"result-gp.dat" u 1:7 t "Yoyo vitesse y" w l,\
"result-gp.dat" u 1:8 t "contrainte" w l


set origin 0.0,0.0
plot\
"result-gp.dat" u 1:9 t "y-h" w l,\
"fichier.dat" u 1:2 t "Controle" w l

set nomultiplot



