function[reached,vL]=robotLimits(v,vmax)
 
 ndof=pa10Ndof();
 vmin=-vmax;
 [vL,zmax] = min(v,vmax);
 [vL,zmin] = max(vL,vmin);

 limits=sum(zmax)+sum(zmin);
 if (limits==2*ndof) then
	reached=%F;
else	
	reached=%T;
end
	
endfunction

//
// Loop
//
function [] =simulationLoop(q0,T,funcControl,funcModel,funcSpy,bSpy,ndof)

global GLOBAL_HIST;

t=0;
// vector x=(q,qdot)
x=[q0,zeros(1,ndof)];


GLOBAL_HIST=list(q0',zeros(ndof,1),zeros(ndof,1),zeros(ndof,1));

while t <= T,

// Robot Control

GLOBAL_U = funcControl(t,x);

// Integration
x=ode(x,t,t+GLOBAL_DT,funcModel);

printf("Time %f\r",t);


// Stockage
GLOBAL_HIST(1)=[GLOBAL_HIST(1),x(1:ndof)']
GLOBAL_HIST(2)=[GLOBAL_HIST(2),x(ndof+1:2*ndof)']
GLOBAL_HIST(3)=[GLOBAL_HIST(3),GLOBAL_U]

// Spy Robot State
if bSpy then
	funcSpy(x);
end;

// Time incrementation,
 t=t+GLOBAL_DT;
 
end; // while

endfunction


