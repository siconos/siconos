// global GLOBAL_U;

//
// Robot Model links
//

function [] = pa10ModelInit()

link('lib/pa10/robotmodel.a',['Ndof','Friction','Inertia','NLEffects', 'Tags'],'c');

endfunction

function [var]=pa10Ndof()
  var=call('Ndof','out',[1,1],1,'i')
endfunction

function [var]=pa10Friction(q,qdot)
  var=call('Friction',q,2,'d',qdot,3,'d','out',[7,1],1,'d')
endfunction

function [var]=pa10Inertia(q)
  var=call('Inertia',q,2,'d','out',[7,7],1,'d')
endfunction

function [var]=pa10NLEffects(q, qdot)
  var=call('NLEffects',q,2,'d',qdot,3,'d','out',[7,1],1,'d')
endfunction

function [var]=pa10Tags(q)
  var=call('Tags',q,2,'d','out',[24,1],1,'d')
endfunction

// Robot Vector Gravity
function [g] = pa10Gravity(q)

dqotz = zeros(q);
[g] = pa10NLEffects(q, dqotz)

endfunction


//
// Robot Lagrange Dynamics
//
function [qddot] = pa10Dynamics(q, qdot, Tau)

// Matrix computation
 M = pa10Inertia(q);

//  motor constants (link-referenced in kg-m^2) 
 M(1,1)=0.75 + M(1,1);
 M(2,2)=0.75 + M(2,2);
 M(3,3)=0.2125 + M(3,3);
 M(4,4)=0.2125 + M(4,4);
 M(5,5)=0.00575 + M(5,5);
 M(6,6)=0.00575 + M(6,6);
 M(7,7)=0.00575 + M(7,7);


 NL = pa10NLEffects(q, qdot);
 
 F = pa10Friction(q,qdot);

// Dynamics equation
[qddot] = inv(M)*(Tau - NL - F);

endfunction

//
// PA10 Model Simulation
//
function [xdot] = pa10Model(t,x)

q=x(1:7)';
qdot=x(8:14)';

qqdot = pa10Dynamics(q, qdot, GLOBAL_U);

xdot=[qdot',qqdot'];

endfunction

//
// PA10 Spying Robot State
//

function [] = pa10Spy(x)

q=x(1:7)';
qmax=[%pi,%pi*97/180,%pi,%pi*143/180,%pi*270/180,%pi,%pi*270/180];

qdot=x(8:14)';
qdotmax=[%pi*57/180,%pi*57/180,%pi*114/180,%pi*114/180,%pi*2,%pi*2,%pi*2]';

tormax=[158,158,68,68,17,17,17]';

ndof=pa10Ndof();

if robotLimits(q,qmax,ndof) then
	  printf("\nerror:: PA10 Joint limit Reached\n");
	  abort;
	end;

if robotLimits(qdot,qdotmax,ndof) then
	printf("\nerror:: PA10 Velocity limit Reached\n");
	abort;
end;

endfunction
