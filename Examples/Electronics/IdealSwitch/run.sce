//initial condition. set alpha
exec("init.sce");

st = 7.45;

exec("computeHandG.sce");

exec("buildMLCP.sce");

// the solution of this first MLCP is :
nextL = [2.999496e-02;
-1.999997e-02;
1.999997e+01;
9.994988e-03;
7.450000e+00;
0.000000e+00;      
0.000000e+00;       
1.999997e+01;       
9.999990e+02];       




exec("computeNext.sce");
exec("computeResidu.sce");
Rr
Ry
Rx


//alpha <-- alpha+1
exec("postStep.sce");

//convergence ok
X_k=X_alpha;
st = 7.4;


exec("computeHandG.sce");

exec("buildMLCP.sce");



nextL=[
3.998494e-02;
-1.999996e-02;
1.999996e+01;
1.998498e-02;
7.400000e+00;
0.000000e+00;
0.000000e+00;
1.999996e+01;
9.999990e+02]



exec("computeNext.sce");

exec("computeResidu.sce");
Rr
Ry
Rx

//alpha <-- alpha+1
exec("postStep.sce");

//convergence ok
X_k=X_alpha;
st = 7.35;


exec("computeHandG.sce");

exec("buildMLCP.sce");

nextL=[4.996992e-02;
-1.999995e-02;
1.999995e+01;
2.996997e-02;
7.350000e+00;
0.000000e+00;
0.000000e+00;
1.999995e+01;
9.999990e+02];

exec("computeNext.sce");

exec("computeResidu.sce");
Rr
Ry
Rx
