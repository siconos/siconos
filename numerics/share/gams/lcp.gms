$if not set filename $set filename 'lcp.gdx'
$if not %gams.user1% == "" $set filename %gams.user1%

$if not set outfile $set outfile 'lcp_sol.gdx'
$if not %gams.user2% == "" $set outfile %gams.user2%

Set i;

parameters M(i, i), q(i), sol(i);

$GDXIN '%filename%';
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

execute_unloadIdx '%outfile%', sol;
