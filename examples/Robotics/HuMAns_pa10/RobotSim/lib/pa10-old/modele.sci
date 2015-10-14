// global GLOBAL_U;

//
// Robot Model links
//

function [] = pa10ModelInit()

link('lib/pa10/robotmodel.a',['modele_coriolis','modele_gravite','modele_inertie','modele_frottements','modele_nddl' ],'c');

endfunction

function [var]=pa10coriolis(q, qdot)
  var=fort('modele_coriolis',q,1,'d',qdot,2,'d','out',[7,7],3,'d')
endfunction

function [var]=pa10Gravity(q)
  var=fort('modele_gravite',q,1,'d','out',[7,1],2,'d')
endfunction

function [var]=pa10inertie(q)
  var=fort('modele_inertie',q,1,'d','out',[7,7],2,'d')
endfunction

function [var]=pa10frottements(q,qdot)
  var=fort('modele_frottements',q,1,'d',qdot,2,'d','out',[7,1],3,'d')
endfunction

function [var]=pa10Ndof()
  var=fort('modele_nddl','out',[1,1],1,'i')
endfunction

function [NL]=arm2NLEffects()
  N = pa10coriolis(q, qdot)*qdot;
  G = pa10Gravity(q);
  [NL]=N+G;
endfunction

//
// Robot Lagrange Dynamics
//
function [qddot] = pa10Dynamics(q, qdot, Tau)

// Matrix computation
 M = pa10inertie(q);
 
 N = pa10coriolis(q, qdot)*qdot;
 G = pa10Gravity(q);
 F = pa10frottements(q,qdot);

// Dynamics equation
//[qddot] = inv(M)*(Tau - N - G - F);
[qddot] = inv(M)*(Tau - G - F);

endfunction

//
// PA10 Model Simulation
//
function [xdot] = pa10Model(t,x)

q=x(1:pa10Ndof())';
qdot=x(pa10Ndof()+1:2*pa10Ndof())';

qqdot = pa10Dynamics(q, qdot, GLOBAL_U);

xdot=[qdot',qqdot'];

endfunction
