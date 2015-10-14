$TEMPLATE[spgvd.real.min_size_iwork.args]
JOBZ,N
$TEMPLATE[spgvd.real.min_size_iwork]
if ( jobz == 'N' || n < 2 )
    return 1;
else
    return 3 + 5*n;
$TEMPLATE[spgvd.real.min_size_work.args]
JOBZ,N
$TEMPLATE[spgvd.real.min_size_work]
if ( n < 2 )
    return 1;
else {
    if ( jobz == 'N' )
        return 2*n;
    else
        return 1 + 6*n + n*n;
}
$TEMPLATE[spgvd.all.UPLO.trait_of]
AP
$TEMPLATE[end]
