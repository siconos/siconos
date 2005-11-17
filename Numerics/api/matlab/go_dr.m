clc;
clear all;

mex -DMEXFLAG -v -g intlab_dr.c ../../src/NSSpack/.libs/libNSSpack.a;

load testlab/Mdr2;
load testlab/qdr2;
qdr2=-qdr2;
load testlab/adr2;
load testlab/bdr2;
bdr2=-bdr2;


n    = int32(40);
itt1  = int32(3000);
itt  = int32(3000);
chat = int32(1);


param1=struct('name','Latin','itermax',itt,'tol',0.000001,'chat',chat,'k_latin',0.007);
param2=struct('name','NLGS','itermax',itt,'tol',0.000001,'chat',chat);


sprintf('\n Latin test \n')


[z1,w1,info1]=intlab_dr(Mdr2,qdr2,n,adr2,bdr2,param1);

sprintf('min(w) %g min(z) %g complementarity %g',min(w1),min(z1),z1'*w1) %'

sprintf(' \n NLGS test \n' )


[z2,w2,info2]=intlab_dr(Mdr2,qdr2,n,adr2,bdr2,param2);

sprintf('min(w) %g min(z) %g complementarity %g',min(w2),min(z2),z2'*w2) %'




