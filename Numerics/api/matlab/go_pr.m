clc;
clear all;

mex -DMEXFLAG -v -g intlab_pr.c ../../src/NSSpack/.libs/libNSSpack.a;

load testlab/Mpr2;
load testlab/qpr2;
qpr2=-qpr2;
load testlab/apr2;
load testlab/bpr2;
bpr2=-bpr2;


n    = int32(40);
itt1  = int32(3000);
itt  = int32(800);
chat = int32(1);


param1=struct('name','Latin','itermax',itt1,'tol',0.000001,'chat',chat,'k_latin',0.007);



param2=struct('name','NLGS','itermax',itt,'tol',0.000001,'chat',chat);


sprintf('\n Latin test \n')


[z1,w1,info1]=intlab_pr(Mpr2,qpr2,n,apr2,bpr2,param1);


sprintf(' \n NLGS test \n' )

[z2,w2,info2]=intlab_pr(Mpr2,qpr2,n,apr2,bpr2,param2);






