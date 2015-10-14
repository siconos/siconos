set label 1 sprintf("Siconos NewtonEulerDS simulation.")  at first 1.75 , first 1.8
#outfile = sprintf('frame%1.jpg',i)
#set output  outfile
plot resultfile1 u 3*i+1:3*i+3 w lp t 'beam1',  resultfile2 u 3*i+1:3*i+3 w lp t 'beam2',  resultfile3 u 3*i+1:3*i+3 w lp t 'beam3'
i=i+every
print "iteration = ", i, " output in animated gif file : ", outputfile
if (i < n) reread