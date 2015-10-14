
sets i, j, p;
parameter M(i,i), H(i,j), E(j,j), Ht(i, j), q(i), b(j), bt(j), Ak(p,j);
$gdxin 'fc3d_avi.gdx'

$load i j p
$load M H E Ht Ak q b bt

$gdxin

display M, H, E, Ht, Ak, q, b, bt

alias(i,k);
alias(j,l);

variables v(i), r(j), y(j);
equations F_v(k), F_r(l), F_y(l), cons_r(p), cons_y(p);

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

reaction(l) = r.l(l)
velocity(l) = v.l(l)

display v.l, r.l, y.l

$exit

* following only works if Q is symmetric
equation defobj; variable obj;
defobj.. obj =e= 0.5*sum(k, x(k)*sum(j, Q(k,j)*x(j))) + sum(k, c(k)*x(k)); 

model qp / defobj, cons /;
solve qp using qcp min obj;
