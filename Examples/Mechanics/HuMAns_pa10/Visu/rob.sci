exec('Load.sci');
m=fscanfMat('../result.dat');
n1=1;
n=size(m,1);
q=zeros(7,n);
q(2:2:7,n1:n)=m(n1:n,2:2:7)';
Visu(q,q)


