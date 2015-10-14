#set  term X11 


toffset = 2.100E-04

basheight = 0.2
heightoff = 0.2
winratio = 1.0
winheight = basheight*winratio

set lmargin 8 
set bmargin 0
set tmargin 0

set size 1.0 , 1.0

set xzeroaxis
set format x ""
 
set grid xtics

set multiplot

set size winratio,winheight

set origin 0.0,winheight*3.0+heightoff
set ylabel "V" 1
plot\
  "SMASHN0p25TRAP1us.dat" u ($1-toffset):2 t "capacitor voltage , SMASH   TRAP 1us" w l,\
  "DiodeBridgeTRAP1us.dat" u 1:2 t           "capacitor voltage , SICONOS TRAP 1us" w l,\
  "SMASHN0p25BE1us.dat" u ($1-toffset):2 t   "capacitor voltage , SMASH   BE 1us" w l,\
  "DiodeBridgeBE1us.dat" u 1:2 t             "capacitor voltage , SICONOS BE 1us" w l
set origin 0.0,winheight*2.0+heightoff
set ylabel "A" 1
plot\
  "SMASHN0p25TRAP1us.dat" u ($1-toffset):3 t "inductor current , SMASH   TRAP 1us" w l,\
  "DiodeBridgeTRAP1us.dat" u 1:3 t           "inductor current , SICONOS TRAP 1us" w l,\
  "SMASHN0p25BE1us.dat" u ($1-toffset):3 t   "inductor current , SMASH   BE 1us" w l,\
  "DiodeBridgeBE1us.dat" u 1:3 t             "inductor current , SICONOS BE 1us" w l
set origin 0.0,winheight+heightoff
set ylabel "V" 1
plot\
  "SMASHN0p25TRAP1us.dat" u ($1-toffset):4 t       "resistor voltage , SMASH   TRAP 1us" w l,\
  "DiodeBridgeTRAP1us.dat" u 1:(-($5 + $6)) t      "resistor voltage , SICONOS TRAP 1us" w l,\
  "SMASHN0p25BE1us.dat" u ($1-toffset):4 t         "resistor voltage , SMASH   BE 1us" w l,\
  "DiodeBridgeBE1us.dat" u 1:(-($5 + $6)) t        "resistor voltage , SICONOS BE 1us" w l

set bmargin 0
set format
set xtics axis
set xlabel "time in s"
set origin 0.0,0.0+heightoff
set ylabel "A" 1
plot\
  "SMASHN0p25TRAP1us.dat" u ($1-toffset):5 t    "resistor current , SMASH   TRAP 1us" w l,\
  "DiodeBridgeTRAP1us.dat" u 1:($4 + $7) t      "resistor current , SICONOS TRAP 1us" w l,\
  "SMASHN0p25BE1us.dat" u ($1-toffset):5 t      "resistor current , SMASH   BE 1us" w l,\
  "DiodeBridgeBE1us.dat" u 1:($4 + $7) t        "resistor current , SICONOS BE 1us" w l
set nomultiplot

