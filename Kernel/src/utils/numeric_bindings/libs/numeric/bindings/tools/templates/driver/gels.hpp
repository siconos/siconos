$TEMPLATE[gels.all.min_size_work.args]
M,N,NRHS
$TEMPLATE[gels.all.min_size_work]
$INTEGER_TYPE minmn = std::min< $INTEGER_TYPE >( m, n );
return std::max< $INTEGER_TYPE >( 1, minmn + std::max< $INTEGER_TYPE >( minmn, nrhs ) );
$TEMPLATE[end]
