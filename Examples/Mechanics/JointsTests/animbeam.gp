#set label 1 sprintf("Siconos ED. gamma = 0.19444 t=%f",i*20*5e-3)  at first 1.75 , first 1.8
#outfile = sprintf('frame%1.jpg',i)
#set output  outfile
plot resultfile1 u 3*i+1:3*i+3 w lp,  resultfile2 u 3*i+1:3*i+3 w lp,  resultfile3 u 3*i+1:3*i+3 w lp 
i=i+every
print i
if (i < n) reread