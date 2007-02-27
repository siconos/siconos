// 
// Pa10 Free in space
//
//
function [U] = pa10Free(t, x)
  U = zeros(pa10Ndof(),1);
endfunction

// 
// Pa10 Joint Control
//
//
function [U] = pa10JointControl(t, x)

global q0;

 q = x(1:pa10Ndof())';
 qdot = x(pa10Ndof()+1:2*pa10Ndof())';

// traj generation
// Constante 
//qd=[0,0.1,0,0,0,0,0]';
// or Sinus
//qd=[0,0.05*sin(0.5*t),0,0 ,0,0,0]';

if (t==0) then
 	q0=q;
end;

 
 qd =f_trajectory(t,q0',[1,1,0,1,0,0,0]',3);

 err = f_error(q, qd);

 verr= f_verror(t,err,0.005,0.6);

 tpd = f_pd(err,verr,[1000,1000,100,1000,50,50,50]',[100,100,10,100,0,0,0]');

 grav =f_gravity(q,pa10Gravity);

 tfric=f_friction(q);

 torque = f_add(tpd, grav,tfric);

 U = f_saturation(torque,[158,158,68,68,17,17,17]');

endfunction


// 
// Pa10 Cartesian Control
//
//
function [U] = pa10CartesianControl(t, x)

global x0;
global xf;

q = x(1:pa10Ndof())';
qdot = x(pa10Ndof()+1:2*pa10Ndof())';
dt = 0.005;

// traj generation
// Constante 
//xd=[0 0 0.1 0 0 0];

if (t==0) then
 	x0=f_xKin(q,pa10Kin,7);
        xf=x0;
   // Test sur X   xf(1)=x0(1)+0.3;
   // Test sur Y   xf(2)=x0(2)+0.3;
   // Test sur Z   
    xf(3)=x0(3)-0.3;
       
end;

xd = f_xtrajectory(t,x0',xf',3);
//xd=x0';

 mKin=f_Kin(q,pa10Kin,pa10Ndof());
 
 mKind=f_Hmat(xd);

 xe=f_xerror(mKin,mKind);

 err = f_Jacinv(q,xe,pa10Jac)

 verr= f_verror(t,err,dt,0.6);

 tpd = f_pd(err,verr,[1000,1000,100,1000,50,50,50]',[100,100,10,100,0,0,0]');

 grav =f_gravity(q,pa10Gravity);

 tfric=f_friction(q);

 torque = f_add(tpd, grav,tfric);

 U = f_saturation(torque,[158,158,68,68,17,17,17]');


endfunction
