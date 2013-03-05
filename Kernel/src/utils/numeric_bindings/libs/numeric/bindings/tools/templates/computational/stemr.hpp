$TEMPLATE[stemr.all.min_size_work.args]
N, JOBZ
$TEMPLATE[stemr.all.min_size_work]
if ( jobz == 'V' ) {
    return std::max< $INTEGER_TYPE >( 1, 18*n );
} else {
    return std::max< $INTEGER_TYPE >( 1, 12*n );
}
$TEMPLATE[stemr.all.min_size_iwork.args]
N,JOBZ
$TEMPLATE[stemr.all.min_size_iwork]
if ( jobz == 'V' ) {
    return std::max< $INTEGER_TYPE >( 1, 10*n );
} else {
    return std::max< $INTEGER_TYPE >( 1, 8*n );
}
$TEMPLATE[end]
