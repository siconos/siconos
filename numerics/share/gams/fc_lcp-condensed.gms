$if not set filename $set filename 'fc_lcp-condensed.gdx'
$if not %gams.user1% == "" $set filename %gams.user1%

$if not set outfile $set outfile 'fc_lcp-condensed_sol.gdx'
$if not %gams.user2% == "" $set outfile %gams.user2%

set j /1 * 2/;

sets i, p;
parameter W(i, i), E(i, i), Wt(i, i), q(i), qt(i), Ak(p, i)

* TODO guess_r(i), guess_y(i), guess_lambda_r(p), guess_lambda_y(p);

$gdxin '%filename%';

$loadIdx W E Wt Ak q qt

* guess_r guess_y guess_lambda_r guess_lambda_y

$gdxin

display Ak, W, E, Wt, q, qt;

display i, p;

alias(i, l);

positive variables l_r(p), l_y(p), s_r(i), s_y(i);
equations  F_r(l), F_y(l), cons_r(p), cons_y(p);

parameter sr(l), sy(l), infos(j);

* cons_r.m(p) = guess_lambda_r(p);
* cons_y.m(p) = guess_lambda_y(p);

F_r(l)..
  sum(i, W(l,i)*s_r(i)) + sum(i, E(l, i)*s_y(i)) - sum(p, Ak(p, l)*l_r(p)) + q(l) =n= 0;

F_y(l)..
  sum(i, Wt(l,i)*s_r(i)) + sum(i, E(l, i)*s_y(i)) - sum(p, Ak(p, l)*l_y(p)) + qt(l) =n= 0;

cons_r(p)..
  sum(i, Ak(p,i)*s_r(i)) =n= 0.;

cons_y(p)..
  sum(i, Ak(p,i)*s_y(i)) =n= 0.;

* parameters r_r(l), r_y(l), r_lr(p), r_ly(p);
* r_r(l) = sum(i, W(l,i)*r.l(i)) + sum(i, E(l, i)*y.l(i)) + q(l) - sum(p, Ak(p,l)*guess_lambda_r(p));
* r_y(l) = sum(i, Wt(l,i)*r.l(i)) + sum(i, E(l, i)*y.l(i)) + qt(l) - sum(p, Ak(p,l)*guess_lambda_y(p));

* r_lr(p) = sum(i, Ak(p,i)*r.l(i));
* r_ly(p) = sum(i, Ak(p,i)*y.l(i));

* display r_r, r_y, r_lr, r_ly;

model lcp / F_r.s_r F_y.s_y cons_r.l_r cons_y.l_y /;

* file fx /"%emp.info%"/;
* putclose fx 'vi F_r r F_y y';

solve lcp using mcp;

sr(i) = s_r.l(i);
sy(i) = s_y.l(i);

infos('1') = lcp.modelstat;
infos('2') = lcp.solvestat;

display lcp.modelstat;
display lcp.solvestat;

execute_unloadIdx '%outfile%', sr, sy, infos;
