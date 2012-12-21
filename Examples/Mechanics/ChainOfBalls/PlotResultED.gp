set autoscale
extension=".tex"

set xlabel 'time (s)' 
set grid 

set term X11 0
set ylabel '$y$' 
plot 	'result.dat' using 1:2  title 'Ball 1' w lines lc rgb "blue",\
     	'result.dat' using 1:3  title 'Ball 2' w lines lc rgb "green",\
	'result.dat' using 1:4  title 'Ball 3' w lines lc rgb "red",\
	'result.dat' using 1:5  title 'Ball 4' w lines lc rgb "black",\
	'result.dat' using 1:6  title 'Ball 5' w lines lc rgb "gold",\
	'result.dat' using 1:7  title 'Ball 6' w lines lc rgb "green",\
	'result.dat' using 1:8  title 'Ball 7' w lines lc rgb "orange",\
	'result.dat' using 1:9  title 'Ball 8' w lines lc rgb "magenta",\
	'result.dat' using 1:10  title 'Ball 9' w lines lc rgb "red",\
	'result.dat' using 1:11  title 'Ball 10' w lines lc rgb "yellow"

set term X11 1
set ylabel '$v$'
plot 	'result.dat' using 1:12  title 'Ball 1' w lines lc rgb "blue",\
     	'result.dat' using 1:13  title 'Ball 2' w lines lc rgb "green",\
	'result.dat' using 1:14  title 'Ball 3' w lines lc rgb "red",\
	'result.dat' using 1:15  title 'Ball 4' w lines lc rgb "black",\
	'result.dat' using 1:16  title 'Ball 5' w lines lc rgb "gold",\
	'result.dat' using 1:17  title 'Ball 6' w lines lc rgb "green",\
	'result.dat' using 1:18  title 'Ball 7' w lines lc rgb "orange",\
	'result.dat' using 1:19  title 'Ball 8' w lines lc rgb "magenta",\
	'result.dat' using 1:20  title 'Ball 9' w lines lc rgb "red",\
	'result.dat' using 1:21  title 'Ball 10' w lines lc rgb "yellow"




