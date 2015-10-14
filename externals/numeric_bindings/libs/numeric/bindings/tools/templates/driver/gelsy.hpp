$TEMPLATE[gelsy.all.min_size_work.args]
M,N,NRHS
$TEMPLATE[gelsy.real.min_size_work]
$INTEGER_TYPE minmn = std::min< $INTEGER_TYPE >( m, n );
return std::max< $INTEGER_TYPE >( 1, std::max< $INTEGER_TYPE >( minmn+3*n+1, 2*minmn+nrhs ));
$TEMPLATE[gelsy.complex.min_size_work]
$INTEGER_TYPE minmn = std::min< $INTEGER_TYPE >( m, n );
return std::max< $INTEGER_TYPE >( 1, std::max< $INTEGER_TYPE >( std::max< $INTEGER_TYPE >( 2*minmn, n+1 ), minmn+nrhs ) );
$TEMPLATE[end]
