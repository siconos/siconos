$TEMPLATE[heevr.complex.min_size_work.args]
N
$TEMPLATE[heevr.complex.min_size_work]
return std::max< $INTEGER_TYPE >( 1, 2*n );
$TEMPLATE[heevr.complex.min_size_rwork.args]
N
$TEMPLATE[heevr.complex.min_size_rwork]
return std::max< $INTEGER_TYPE >( 1, 24*n );
$TEMPLATE[heevr.complex.min_size_iwork.args]
N
$TEMPLATE[heevr.complex.min_size_iwork]
return std::max< $INTEGER_TYPE >( 1, 10*n );
$TEMPLATE[end]
