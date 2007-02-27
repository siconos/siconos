// See doc. arm2Control.pdf

function [U] = arm2JointControl(t, x)


q = x(1:2)';

qdot = x(3:4)';

qd = f_trajectory(t);

err = f_error(q, qd);

verr= f_verror(t,err,0.1,0.6);

tpd = f_pd(err,verr,[300.0,500.0]',[100.0,20.0]');

//grav=arm2gravite(q);
grav =f_gravity(q,arm2Gravity);


tfric=f_friction(q);

torque = f_add(tpd, grav,tfric);

U = f_saturation(torque,[800,800]');

//U=[0,0]';

endfunction

