$TEMPLATE[hbevd.complex.min_size_iwork.args]
JOBZ,N
$TEMPLATE[hbevd.complex.min_size_iwork]
if ( jobz == 'N' || n < 2 )
    return 1;
else
    return 3 + 5*n;
$TEMPLATE[hbevd.complex.min_size_rwork.args]
JOBZ,N
$TEMPLATE[hbevd.complex.min_size_rwork]
if ( n < 2 )
    return 1;
else {
    if ( jobz == 'N' )
        return n;
    else
        return 1 + 5*n + 2*n*n;
}
$TEMPLATE[hbevd.complex.min_size_work.args]
JOBZ,N
$TEMPLATE[hbevd.complex.min_size_work]
if ( n < 2 )
    return 1;
else {
    if ( jobz == 'N' )
        return n;
    else
        return 2*n*n;
}
$TEMPLATE[hbevd.all.UPLO.trait_of]
AB
$TEMPLATE[end]
