$TEMPLATE[stevd.real.min_size_iwork.args]
JOBZ,N
$TEMPLATE[stevd.real.min_size_iwork]
if ( jobz == 'N' || n < 2 )
    return 1;
else
    return 3 + 5*n;
$TEMPLATE[stevd.real.min_size_work.args]
JOBZ,N
$TEMPLATE[stevd.real.min_size_work]
if ( jobz == 'N' || n < 2 )
    return 1;
else
    return 1 + 4*n + n*n;
$TEMPLATE[end]
