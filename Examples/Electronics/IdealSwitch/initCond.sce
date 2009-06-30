R1=1e-3;
R2=1e3;
h=1e-7;
sL = 200e-6;


X_k=0.0;
X_alpha=X_k;
L_alpha=[20.0/(R1+R2);
-20.0/(R1+R2);
20.0-20.0*(R1/(R1+R2));
0.0;
7.5;
0.0;
0.0;
20.0-20.0*(R1/(R1+R2));
R2-R1];




st=7.5;


h_alpha=[L_alpha(5)-st;
X_alpha-L_alpha(4);
L_alpha(3)-20+L_alpha(1)*(L_alpha(7)+R1);
L_alpha(3)+L_alpha(2)*(L_alpha(9)+R1);
X_alpha-L_alpha(1)-L_alpha(2);
R2-L_alpha(7)-R1;
L_alpha(5)-L_alpha(4)+L_alpha(6);
R2-L_alpha(9)-R1;
-L_alpha(3)+L_alpha(8)];

g_alpha = (L_alpha(3)-L_alpha(4))/sL;