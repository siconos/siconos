set autoscale
extension=".tex"
resultfile='result.dat'
#resultfile="resultMonodisperseChain_LZBModel.ref"
#resultfile="resultTaperedChain_LZBModel.ref"
apple=1 



set xlabel 'time (s)' 
set grid 


if (apple == 1) {
set term aqua 0;
} else {
set term x11 0;
}

set ylabel '$y$' 
plot 	resultfile every::2 using 1:2  title 'Ball 1' w lines lc rgb "blue",\
     	resultfile every::2 using 1:3  title 'Ball 2' w lines lc rgb "green",\
	resultfile every::2 using 1:4  title 'Ball 3' w lines lc rgb "red",\
	resultfile every::2 using 1:5  title 'Ball 4' w lines lc rgb "black",\
	resultfile every::2 using 1:6  title 'Ball 5' w lines lc rgb "gold",\
	resultfile every::2 using 1:7  title 'Ball 6' w lines lc rgb "green",\
	resultfile every::2 using 1:8  title 'Ball 7' w lines lc rgb "orange",\
	resultfile every::2 using 1:9  title 'Ball 8' w lines lc rgb "magenta",\
	resultfile every::2 using 1:10  title 'Ball 9' w lines lc rgb "red",\
	resultfile every::2 using 1:11  title 'Ball 10' w lines lc rgb "yellow"

if (apple == 1) {
set term aqua 1
} else {
set term x11 1;
}

set ylabel '$v$'
plot 	resultfile every::2 using 1:12  title 'Ball 1' w lines lc rgb "blue",\
     	resultfile every::2 using 1:13  title 'Ball 2' w lines lc rgb "green",\
	resultfile every::2 using 1:14  title 'Ball 3' w lines lc rgb "red",\
	resultfile every::2 using 1:15  title 'Ball 4' w lines lc rgb "black",\
	resultfile every::2 using 1:16  title 'Ball 5' w lines lc rgb "gold",\
	resultfile every::2 using 1:17  title 'Ball 6' w lines lc rgb "green",\
	resultfile every::2 using 1:18  title 'Ball 7' w lines lc rgb "orange",\
	resultfile every::2 using 1:19  title 'Ball 8' w lines lc rgb "magenta",\
	resultfile every::2 using 1:20  title 'Ball 9' w lines lc rgb "red",\
	resultfile every::2 using 1:21  title 'Ball 10' w lines lc rgb "yellow"




