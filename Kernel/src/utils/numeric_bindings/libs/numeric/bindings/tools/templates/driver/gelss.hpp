$TEMPLATE[gelss.real.min_size_work.args]
M,N,NRHS
$TEMPLATE[gelss.real.min_size_work]
$INTEGER_TYPE minmn = std::min< $INTEGER_TYPE >( m, n );
return std::max< $INTEGER_TYPE >( 1, 3*minmn + std::max< $INTEGER_TYPE >( std::max< $INTEGER_TYPE >( 2*minmn, std::max< $INTEGER_TYPE >(m,n) ), nrhs ) );
$TEMPLATE[gelss.complex.extra_variables]
MINMN
$TEMPLATE[gelss.complex.extra_opt_variables]
MINMN
$TEMPLATE[gelss.complex.MINMN.init]
$INTEGER_TYPE minmn = std::min< $INTEGER_TYPE >( size_row(a), size_column(a) );
$TEMPLATE[gelss.complex.min_size_work.args]
M,N,NRHS,MINMN
$TEMPLATE[gelss.complex.min_size_work]
return std::max< $INTEGER_TYPE >( 1, 2*minmn + std::max< $INTEGER_TYPE >( std::max< $INTEGER_TYPE >( m,n ), nrhs ) );
$TEMPLATE[gelss.complex.min_size_rwork.args]
MINMN
$TEMPLATE[gelss.complex.min_size_rwork]
return 5*minmn;
$TEMPLATE[end]
