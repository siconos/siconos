$if not set filename $set filename 'fc3d_avi-condensed.gdx'
$if not %gams.user1% == "" $set filename %gams.user1%

$if not set outfile $set outfile 'fc3d_avi-condensed_sol.gdx'
$if not %gams.user2% == "" $set outfile %gams.user2%

$if not set additional_constr $set additional_constr 'none'
$if not "%gams.user3%" == "" $set additional_constr %gams.user3%

set j /1 * 4/;

sets i, p;
parameter W(i,i), E(i,i), Wt(i, i), q(i), qt(i), Ak(p,i), guess_r(i), guess_y(i), guess_lambda_r(p), guess_lambda_y(p);

$gdxin '%filename%';

$loadIdx W E Wt Ak q qt guess_r guess_y guess_lambda_r guess_lambda_y

$gdxin

*display Ak, W, E, Wt, q, qt;

*display i, p;

alias(i,l);

variables r(l), y(l);

$ifthen %additional_constr% == 'normal_pos'
  r.lo(l)$(mod(ord(l), 3) = 1) = 0.;
  y.lo(l)$(mod(ord(l), 3) = 1) = 0.;
$endif

equations  F_r(l), F_y(l), cons_r(p), cons_y(p);

parameter reaction(l), velocity(l), infos(j);

r.l(i) = guess_r(i);
y.l(i) = guess_y(i);

F_r(l)..
  sum(i, W(l,i)*r(i)) + sum(i, E(l, i)*y(i)) + q(l) =n= 0;

F_y(l)..
  sum(i, Wt(l,i)*r(i)) + sum(i, E(l, i)*y(i)) + qt(l) =n= 0;

cons_r(p)..
  sum(i, Ak(p,i)*r(i)) =g= 0.;

cons_y(p)..
  sum(i, Ak(p,i)*y(i)) =g= 0.;

cons_r.m(p) = guess_lambda_r(p);
cons_y.m(p) = guess_lambda_y(p);

parameters r_r(l), r_y(l), r_lr(p), r_ly(p);
r_r(l) = sum(i, W(l,i)*r.l(i)) + sum(i, E(l, i)*y.l(i)) + q(l) - sum(p, Ak(p,l)*guess_lambda_r(p));
r_y(l) = sum(i, Wt(l,i)*r.l(i)) + sum(i, E(l, i)*y.l(i)) + qt(l) - sum(p, Ak(p,l)*guess_lambda_y(p));

r_lr(p) = sum(i, Ak(p,i)*r.l(i));
r_ly(p) = sum(i, Ak(p,i)*y.l(i));

*display r_r, r_y, r_lr, r_ly;

model vi / all /;

file fx /"%emp.info%"/;
putclose fx 'vi F_r r F_y y';

solve vi using emp;

reaction(i) = r.l(i);
velocity(i) = sum(l, W(i,l)*r.l(l)) + q(i) + sum(l, E(i, l)*y.l(l));

infos('1') = vi.modelstat;
infos('2') = vi.solvestat;
infos('3') = vi.iterusd;
infos('4') = vi.Resusd;

execute_unloadIdx '%outfile%', reaction, velocity, infos;
