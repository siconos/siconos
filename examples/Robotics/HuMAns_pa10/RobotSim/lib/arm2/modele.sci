// global GLOBAL_U;

//
// Robot Model links
//

function [] = arm2ModelInit()

link('lib/arm2/robotmodel.a',['Ndof','Friction','Inertia','NLEffects'],'c');

endfunction

function [var]=arm2Ndof()
  var=fort('Ndof','out',[1,1],1,'i')
endfunction

function [var]=arm2Friction(q,qdot)
  var=fort('Friction',q,2,'d',qdot,3,'d','out',[2,1],1,'d')
endfunction

function [var]=arm2Inertia(q)
  var=fort('Inertia',q,2,'d','out',[2,2],1,'d')
endfunction

function [var]=arm2NLEffects(q, qdot)
  var=call('NLEffects',q,2,'d',qdot,3,'d','out',[2,1],1,'d')
endfunction

function [var]=arm2Tags(q)
  var=call('Tags',q,2,'d','out',[6,1],1,'d')
endfunction

// Robot Vector Gravity
function [g] = arm2Gravity(q)

dqotz = zeros(q);
[g] = arm2NLEffects(q, dqotz)

endfunction


//
// Robot Lagrange Dynamics
//
function [qddot] = arm2Dynamics(q, qdot, Tau)

// Matrix computation
 M = arm2Inertia(q);
 
 NL = arm2NLEffects(q, qdot);
 
 F = arm2Friction(q,qdot);

// Dynamics equation
[qddot] = inv(M)*(Tau - NL - F);

endfunction

//
// Arm2 Model Simulation
//
function [xdot] = arm2Model(t,x)

q=x(1:2)';
qdot=x(3:4)';

qqdot = arm2Dynamics(q, qdot, GLOBAL_U);

xdot=[qdot',qqdot'];

endfunction
