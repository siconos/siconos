$if not set filename $set filename 'fc3d_avi.gdx'
$if not "%gams.user1%" == "" $set filename %gams.user1%

$if not set outfile $set outfile 'fc3d_avi_sol.gdx'
$if not "%gams.user2%" == "" $set outfile %gams.user2%

$if not set additional_constr $set additional_constr 'none'
$if not "%gams.user3%" == "" $set additional_constr %gams.user3%


set info_l /1 * 4/;

sets i, j, p;
parameter M(i,i), H(i,j), E(j,j), Ht(i, j), q(i), b(j), bt(j), Ak(p,j);
$gdxin '%filename%'

$loadIdx M H E Ht Ak q b bt

$gdxin


alias(i,k);
alias(j,l);

variables v(i), r(j), y(j);
equations F_v(k), F_r(l), F_y(l), cons_r(p), cons_y(p);

parameter reaction(l), velocity(i), infos(info_l);

$ifthen %additional_constr% == 'normal_pos'
  r.lo(l)$(mod(ord(l), 3) = 1) = 0.;
  y.lo(l)$(mod(ord(l), 3) = 1) = 0.;
$endif

F_v(k)..
  sum(i, M(k,i)*v(i)) + sum(j, -H(k,j)*r(j)) - q(k) =n= 0;

F_r(l)..
  sum(i, H(i,l)*v(i)) + sum(j, E(l, j)*y(j)) + b(l) =n= 0;

F_y(l)..
  sum(i, Ht(i,l)*v(i)) + sum(j, E(l, j)*y(j)) + bt(l) =n= 0;

cons_r(p)..
  sum(j, Ak(p,j)*r(j)) =g= 0.;

cons_y(p)..
  sum(j, Ak(p,j)*y(j)) =g= 0.;

model vi / all /;

file fx /"%emp.info%"/;
putclose fx 'vi F_v v F_r r F_y y';

* x.up(j) = 1e2;

solve vi using emp;

reaction(l) = r.l(l);
velocity(i) = v.l(i);

display v.l, r.l, y.l;

infos('1') = vi.modelstat;
infos('2') = vi.solvestat;
infos('3') = vi.iterusd;
infos('4') = vi.Resusd;

execute_unloadIdx '%outfile%', reaction, velocity, infos;
