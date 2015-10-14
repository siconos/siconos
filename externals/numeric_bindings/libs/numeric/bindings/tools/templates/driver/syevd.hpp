$TEMPLATE[syevd.real.min_size_iwork.args]
JOBZ,N
$TEMPLATE[syevd.real.min_size_iwork]
if ( jobz == 'N' || n < 2 )
    return 1;
else
    return 3 + 5*n;
$TEMPLATE[syevd.real.min_size_work.args]
JOBZ,N
$TEMPLATE[syevd.real.min_size_work]
if ( n < 2 )
    return 1;
else {
    if ( jobz == 'N' )
        return 2*n + 1;
    else
        return 1 + 6*n + 2*n*n;
}
$TEMPLATE[end]
