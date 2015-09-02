
set j /1 * 2/;

sets i, p;
parameter W(i,i), E(i,i), Wt(i, i), q(i), qt(i), Ak(p,i);

$gdxin 'fc3d_avi-condensed.gdx'

$loadIdx W E Wt Ak q qt

$gdxin

display Ak, W, E, Wt, q, qt;

alias(i,l);

variables r(l), y(l);
equations  F_r(l), F_y(l), cons_r(p), cons_y(p);

parameter reaction(l), velocity(l), infos(j);

F_r(l)..
  sum(i, W(l,i)*r(i)) + sum(i, E(l, i)*y(i)) + q(l) =n= 0;

F_y(l)..
  sum(i, Wt(l,i)*r(i)) + sum(i, E(l, i)*y(i)) + qt(l) =n= 0;

cons_r(p)..
  sum(i, Ak(p,i)*r(i)) =g= 0.;

cons_y(p)..
  sum(i, Ak(p,i)*y(i)) =g= 0.;

model vi / all /;

file fx /"%emp.info%"/;
putclose fx 'vi F_r r F_y y';

solve vi using emp;

reaction(i) = r.l(i);
velocity(i) = r.m(i);

* display r.l, y.l

infos('1') = vi.modelstat;
infos('2') = vi.solvestat;

display vi.modelstat;
display vi.solvestat;

execute_unloadIdx 'fc3d_avi-condensed_sol.gdx', reaction, velocity, infos;
