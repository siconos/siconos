$TEMPLATE[steqr.all.N.trait]
size,D
$TEMPLATE[steqr.all.min_size_work.args]
N, COMPZ
$TEMPLATE[steqr.all.min_size_work]
if ( compz == 'N' ) {
    return 1;
} else {
    return std::max< $INTEGER_TYPE >( 1, 2*n-2 );
}
$TEMPLATE[end]
