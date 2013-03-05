$TEMPLATE[sbgvd.real.min_size_iwork.args]
JOBZ,N
$TEMPLATE[sbgvd.real.min_size_iwork]
if ( jobz == 'N' || n < 2 )
    return 1;
else
    return 3 + 5*n;
$TEMPLATE[sbgvd.real.min_size_work.args]
JOBZ,N
$TEMPLATE[sbgvd.real.min_size_work]
if ( n < 2 )
    return 1;
else {
    if ( jobz == 'N' )
        return 3*n;
    else
        return 1 + 5*n + 2*n*n;
}
$TEMPLATE[sbgvd.all.UPLO.trait_of]
AB
$TEMPLATE[end]
