
clear all;

mex -DMEXFLAG -v -g intlab_lcp.c  ../../src/NSSpack/.libs/libNSSpack.a;

load testlab/Mpoin;
load testlab/qpoin;


n    = int32(31);
itt  = int32(2000);
chat = int32(1);


param1=struct('name','Latin','itermax',itt,'tol',0.00000000001,'chat',chat,'k_latin',0.007);

param2=struct('name','NLGS','itermax',itt,'tol',0.00000000001,'chat',chat);

param3=struct('name','CPG','itermax',itt,'tol',0.00000000001,'chat',chat);

param4=struct('name','LexicoLemke','itermax',itt,'chat',chat);

param5=struct('name','QP','tol',0.00000000001,'chat',chat);
param6=struct('name','NSQP','tol',0.00000000001,'chat',chat);

param7=struct('name','Latin_w','itermax',itt,'tol',0.00000000001,'chat',chat,'k_latin',0.007,'relax',0.5);

param8=struct('name','NewtonMin','itermax',itt,'tol',0.00000000001,'chat',chat);







sprintf('\n Latin test \n')


[z1,w1,info1]=intlab_lcp(Mpoin,qpoin,n,param1);

sprintf('min(w) %g min(z) %g complementarity %g',min(w1),min(z1),z1'*w1) %'

sprintf(' \n NLGS test \n' )


[z2,w2,info2]=intlab_lcp(Mpoin,qpoin,n,param2);

sprintf('min(w) %g min(z) %g complementarity %g',min(w2),min(z2),z2'*w2) %'


sprintf('\n CPG test \n')


[z3,w3,info3]=intlab_lcp(Mpoin,qpoin,n,param3);

sprintf('min(w) %g min(z) %g complementarity %g',min(w3),min(z3),z3'*w3) %'


sprintf('\n Lexicolemke test \n')


[z4,w4,info4]=intlab_lcp(Mpoin,qpoin,n,param4);

sprintf('min(w) %g min(z) %g complementarity %g',min(w4),min(z4),z4'*w4) %'

sprintf('\n QP test \n')

[z5,w5,info5]=intlab_lcp(Mpoin,qpoin,n,param5);

sprintf('min(w) %g min(z) %g complementarity %g',min(w5),min(z5),z5'*w5) %'

sprintf('\n NSQP test \n')

    
[z6,w6,info6]=intlab_lcp(Mpoin,qpoin,n,param6);

sprintf('min(w) %g min(z) %g complementarity %g',min(w6),min(z6),z6'*w6) %'

sprintf('\n Latin_w test \n')


[z7,w7,info7]=intlab_lcp(Mpoin,qpoin,n,param7);

sprintf('min(w) %g min(z) %g complementarity %g',min(w7),min(z7),z7'*w7) %'


sprintf('\n NewtonMin test \n')


[z8,w8,info8]=intlab_lcp(Mpoin,qpoin,n,param8);

sprintf('min(w) %g min(z) %g complementarity %g',min(w8),min(z8),z8'*w8) %'






