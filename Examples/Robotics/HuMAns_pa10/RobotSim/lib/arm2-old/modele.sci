// global GLOBAL_U;

//
// Robot Model links
//

function [] = arm2ModelInit()

link('lib/arm2/robotmodel.a',['modele_coriolis','modele_gravite','modele_inertie','modele_frottements','modele_nddl' ],'c');

endfunction

function [var]=arm2coriolis(q, qdot)
  var=call('modele_coriolis',q,1,'d',qdot,2,'d','out',[2,2],3,'d')
endfunction

function [var]=arm2Gravity(q)
  var=fort('modele_gravite',q,1,'d','out',[2,1],2,'d')
endfunction

function [var]=arm2inertie(q)
  var=fort('modele_inertie',q,1,'d','out',[2,2],2,'d')
endfunction

function [var]=arm2frottements(q,qdot)
  var=fort('modele_frottements',q,1,'d',qdot,2,'d','out',[2,1],3,'d')
endfunction

function [var]=arm2Ndof()
  var=fort('modele_nddl','out',[1,1],1,'i')
endfunction

function [NL]=arm2NLEffects()
  N = arm2coriolis(q, qdot)*qdot;
  G = arm2Gravity(q);
  [NL]=N+G;
endfunction

//
// Robot Lagrange Dynamics
//
function [qddot] = arm2Dynamics(q, qdot, Tau)

// Matrix computation
 M = arm2inertie(q);
 
 N = arm2coriolis(q, qdot)*qdot;
 G = arm2Gravity(q);
 F = arm2frottements(q,qdot);

// Dynamics equation
[qddot] = inv(M)*(Tau - N - G - F);

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
