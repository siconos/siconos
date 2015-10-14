$TEMPLATE[bdsdc.real.min_size_work.args]
COMPQ, N
$TEMPLATE[bdsdc.real.min_size_work]
switch ( compq ) {
    case 'N': return 4*n;
    case 'P': return 6*n;
    case 'I': return 3*n*n + 4*n;
}
$TEMPLATE[end]
