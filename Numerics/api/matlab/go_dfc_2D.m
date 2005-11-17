clc;
clear all;



mex  -DMEXFLAG -v -g intlab_dfc_2D.c  ../../src/NSSpack/.libs/libNSSpack.a;




load testlab/K1_2;


vec=zeros(796*796,1);

for i=1:796
for j =1:796
vec(i+(j-1)*796)=K1(i,j);
end
end


%keyboard
load testlab/F1_2;
load testlab/J1_2;

load testlab/ddl_n_2;
load testlab/ddl_t_2;
load testlab/ddl_d_2;

ddl_n=int32(ddl_n');%'
ddl_t=int32(ddl_t');%'
ddl_d=int32(ddl_d');%'




n      = int32(796);
dim_t  = int32(25);
dim_d  = int32(17);
itt    = int32(1000);
chat   = int32(1);
mu     = 0.5;

param1=struct('name','Cfd_latin','itermax',itt,'tol',0.000001,'k_latin',0.6,'chat',chat,'J1',J1,'dim_t',dim_t,'ddl_n',ddl_n,'ddl_t',ddl_t,'dim_d',dim_d,'ddl_d',ddl_d);

param2=struct('name','NSQP','tol',0.000001,'chat',chat,'J1',J1,'dim_t',dim_t,'ddl_n',ddl_n,'ddl_t',ddl_t,'dim_d',dim_d,'ddl_d',ddl_d);

param3=struct('name','Lemke','itermax',itt,'chat',chat,'J1',J1,'dim_t',dim_t,'ddl_n',ddl_n,'ddl_t',ddl_t,'dim_d',dim_d,'ddl_d',ddl_d);

param4=struct('name','NLGS','itermax',itt,'tol',0.000001,'chat',chat,'J1',J1,'dim_t',dim_t,'ddl_n',ddl_n,'ddl_t',ddl_t,'dim_d',dim_d,'ddl_d',ddl_d);


sprintf('\n Cfd_latin test \n')

[z1,w1,info1]=intlab_dfc_2D(vec,F1,n,mu,param1);


sprintf(' \n NSQP test \n' )


[z2,w2,info2]=intlab_dfc_2D(vec,F1,n,mu,param2);


sprintf('\n Lemke test \n')


[z3,w3,info3]=intlab_dfc_2D(vec,F1,n,mu,param3);

sprintf('\n NLGS test \n')


[z4,w4,info4]=intlab_dfc_2D(vec,F1,n,mu,param4);



