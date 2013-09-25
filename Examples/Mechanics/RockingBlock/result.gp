#set  term X11
plot 'result.dat'  using 1:2  t 'x' w lines,\
     'result.dat'  using 1:3  t 'y' w lines,\
     'result.dat'  using 1:4  t 'theta' w lines