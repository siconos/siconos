
Set i;

parameters M(i, i), q(i), sol(i);

$GDXIN lcp.gdx
$loadIdx M, q
$GDXIN

alias(i,j);

positive variable x(i);
equation f(i);

f(i)..
    sum(j, M(i,j)*x(j)) + q(i) =g= 0;

model foo /f.x/;

x.l(i) = 0.1;

solve foo using mcp;

sol(i) = x.l(i)

execute_unloadIdx 'lcp_sol.gdx', sol;
