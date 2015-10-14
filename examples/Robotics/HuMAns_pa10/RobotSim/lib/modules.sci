
function [qd] = f_trajectory_lin(t,qi,qf,duration)

if t < duration then 
    qd = qi + ((t/duration) * (qf - qi));
else
    qd = qf;
end

function [qd] = f_trajectory(t,qi,qf,duration)

if (t==0) then
 d2=duration*duration;
 a2=(3/(d2))*(qf-qi);
 a3=(-2/(d2*duration))*(qf-qi);
end;

if t < duration then 
    t2=t*t;
    qd = qi + a2*t2 + a3*t*t2;
else
    qd = qf;
end
endfunction

function [err] = f_error(q,qd)
 global GLOBAL_HIST;

 err = qd - q;

 GLOBAL_HIST(4)=[GLOBAL_HIST(4),qd];

endfunction

function [verr] = f_verror(t,err,pDt,pAlpha)

global errOld;
global verrOld;

// Init
 if (t==0) then
   errOld = err;
   verrOld= zeros(err);
 end

 verr=(err - errOld)/pDt;

 verr=verr*pAlpha + (1-pAlpha)*verrOld;

errOld=err;
verrOld=verr;

endfunction

function [grav] = f_gravity(q,funcGravity)
 grav = funcGravity(q);
endfunction

function [tpd] = f_pd(err,verr,pKp,pKd)
for i = 1:size(err,1)
 tpd(i) = pKp(i) * err(i) + pKd(i)*verr(i);
end;
endfunction

function [somme] = f_add(input1,input2,input3)
 somme=input1+input2+input3;
endfunction


function [vOutput] = f_saturation(vInput,pSat)
	[reach,vOutput]=robotLimits(vInput,pSat);
//	if reach then
//.		printf("\nWarning:: torque saturation\n");
//	end;
endfunction

function [tfric] = f_friction(q)
// No friction compensation
 tfric=zeros(q);
endfunction

//
// qdot = pseudoinverse(J) xdot
//
function [qdot] = f_Jacinv(q,xdot,funcJac)

pseudoJ=pinv(funcJac(q));

if rank(pseudoJ)<>6 then
	  printf("\nerror:: Singularities Reached\n");
	  abort;
	end;

qdot=pseudoJ*xdot;

endfunction

//
// 
//
function [mkin] = f_Kin(q,funcGeom,index)

 model=funcGeom(q);
 mkin=model(index);

endfunction

//
// 
//
function [x] = f_xKin(q,funcGeom,index)

// Robot Kin
 model=funcGeom(q);
 mkin=model(index);

 x=[0 0 0 0 0 0];
 
 // Translation
 x(1)=mkin(1,4);
 x(2)=mkin(2,4);
 x(3)=mkin(3,4);
 // Euler Rotation
 [x(4),x(5),x(6)]=euler_rot(mkin);


endfunction

//
// 
//
function [mkin] = f_Hmat(x)

 mkin=[0 0 0 0; 0 0 0 0; 0 0 0 0;0 0 0 1];
 // Translation
 mkin(1,4)=x(1);
 mkin(2,4)=x(2);
 mkin(3,4)=x(3);
 // Rotation
 mkin([1 2 3],[1 2 3])=rot_euler(x(4),x(5),x(6));

endfunction

//
//
//
function [qd] = f_xtrajectory(t,qi,qf,duration)


global a2;
global a3;
global d2;

if (t==0) then
 d2=duration*duration;
 a2=(3/(d2))*(qf-qi);
 a3=(-2/(d2*duration))*(qf-qi);
end;

if t < duration then 
    t2=t*t;
    qd = qi + a2*t2 + a3*t*t2;
else
    qd = qf;
end

endfunction


//
// 
//
function [err] = f_xerror(mKin,mKind)

 global GLOBAL_HIST;

err=zeros(6,1);	

// Error in translation
err([1 2 3])=mKind([1 2 3],4) - mKin([1 2 3],4);

// Error in rotation
rot=mKind([1 2 3],[1 2 3])*mKin([1 2 3],[1 2 3])';
err(4) = 0.5*( rot(3,2) - rot(2,3));
err(5) = 0.5*( rot(1,3) - rot(3,1));
err(6) = 0.5*( rot(2,1) - rot(1,2));

 GLOBAL_HIST(4)=[GLOBAL_HIST(4),[err;0]];

endfunction

// compute the rotation matrix 
// using yaw pitch roll parametrization
// Y Roulis : rot(z) phi
// P Tangage: rot(y) theta
// L Lacet  : rot(x) psi
//

function [MROT] = rot_ypl(phi,theta,psi)

cphi=cos(phi);
cthe=cos(theta);
cpsi=cos(psi);
sphi=sin(phi);
sthe=sin(theta);
spsi=sin(psi);

MROT = [[ cphi*cthe cphi*sthe*spsi - sphi*cpsi cphi*sthe*cpsi+sphi*spsi];
        [ sphi*cthe sphi*sthe*spsi + cphi*cpsi sphi*sthe*cpsi-cphi*spsi];
        [ -sthe cthe*spsi cthe*cpsi]
       ];

endfunction;

// compute yaw pitch roll parametrization
// using the rotation matrix
// Y Roulis : rot(z) phi
// P Tangage: rot(y) theta
// L Lacet  : rot(x) psi
//

function [Y,P,L] = ypl_rot(mrot)

sx=mrot(1,1);
nx=mrot(1,2);
ax=mrot(1,3);
sy=mrot(2,1);
ny=mrot(2,2);
ay=mrot(2,3);
sz=mrot(3,1);
nz=mrot(3,2);
az=mrot(3,3);

if ( (sy==0) & (sx==0) ) then
	
end;

Y = atan(sy,sx);
cphi=cos(Y);
sphi=sin(Y);

P = atan(-sz,cphi*sx+sphi*sy);
L = atan(sphi*ax-cphi*ay,-sphi*nx+cphi*ny);

endfunction;

// compute the rotation matrix 
// using euler angle parametrization
// rot(x) phi
// rot(y) theta
// rot(z) psi
//

function [MROT] = rot_euler(phi,theta,psi)

cphi=cos(phi);
cthe=cos(theta);
cpsi=cos(psi);
sphi=sin(phi);
sthe=sin(theta);
spsi=sin(psi);

MROT = [[ cthe*cpsi  cthe*spsi -sthe];
        [ sphi*sthe*cpsi - cphi*spsi  sphi*sthe*spsi+cphi*cpsi sphi*cthe];
        [ cphi*sthe*cpsi+sthe*spsi cphi*sthe*spsi  - sphi*cpsi cphi*cthe];
       ];

endfunction;


// compute euler angle parametrization 
// using the rotation matrix
// rot(z0) phi
// rot(x) theta
// rot(z) psi
//

function [phi,theta,psi] = euler_rot(mrot)

 phi   = atan(mrot(2,3),mrot(3,3));
 theta = asin(-mrot(1,3));
 psi   = atan(mrot(1,2),mrot(1,1));

endfunction;

//           
// Calcule qd_k+1=qk +pDt*pseudoInv(J(q_k))xdotd;
//  
//
function [qk1] = f_Vel(t,pDt,qk,xdot,funcJac)

global GLOBAL_HIST;
global q_old;

// Init
 if (t==0) then
  q_old=qk;
 end

 pseudoJ=pinv(funcJac(qk));

if rank(pseudoJ)<>6 then
	  printf("\nerror:: Singularities Reached\n");
	  abort;
	end;

 qdot=pseudoJ*xdot';

 qk1=qdot*pDt+q_old;

 q_old=qk1;

 GLOBAL_HIST(4)=[GLOBAL_HIST(4),qk1];

endfunction
