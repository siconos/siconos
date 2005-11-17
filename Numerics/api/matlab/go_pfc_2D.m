clc;
clear all;

mex -DMEXFLAG -v -g intlab_pfc_2D.c  ../../src/NSSpack/.libs/libNSSpack.a;

load testlab/W20mu;
load testlab/b20mu;

b20mu=-b20mu';%'

n    = int32(78);
itt  = int32(1000);
chat = int32(1);
mu   = 0.3;


param1=struct('name','Latin','itermax',itt,'tol',0.00000000001,'chat',chat,'k_latin',35.0);
param2=struct('name', 'NLGS','itermax',itt,'tol',0.00000000001,'chat',chat);
param3=struct('name' , 'CPG','itermax',itt,'tol',0.00000000001,'chat',chat);


sprintf('\n Latin test \n')


[z1,w1,info1]=intlab_pfc_2D(W20mu,b20mu,n,mu,param1);


sprintf(' \n NLGS test \n' )


[z2,w2,info2]=intlab_pfc_2D(W20mu,b20mu,n,mu,param2);


sprintf('\n CPG test \n')


[z3,w3,info3]=intlab_pfc_2D(W20mu,b20mu,n,mu,param3);



