$TEMPLATE[sbevd.real.min_size_iwork.args]
JOBZ,N
$TEMPLATE[sbevd.real.min_size_iwork]
if ( jobz == 'N' || n < 2 )
    return 1;
else
    return 3 + 5*n;
$TEMPLATE[sbevd.real.min_size_work.args]
JOBZ,N
$TEMPLATE[sbevd.real.min_size_work]
if ( n < 2 )
    return 1;
else {
    if ( jobz == 'N' )
        return 2*n;
    else
        return 1 + 5*n + 2*n*n;
}
$TEMPLATE[sbevd.all.UPLO.trait_of]
AB
$TEMPLATE[sbevd.all.LIWORK.trait]
size,IWORK
$TEMPLATE[end]
