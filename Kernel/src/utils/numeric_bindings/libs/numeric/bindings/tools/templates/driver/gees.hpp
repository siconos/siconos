$TEMPLATE[gees.real.min_size_work.args]
N
$TEMPLATE[gees.real.min_size_work]
return std::max< $INTEGER_TYPE >( 1, 3*n );
$TEMPLATE[gees.all.min_size_bwork.args]
N, SORT
$TEMPLATE[gees.all.min_size_bwork]
if ( sort == 'N' )
    return 0;
else
    return n;
$TEMPLATE[gees.complex.min_size_work.args]
N
$TEMPLATE[gees.complex.min_size_work]
return std::max< $INTEGER_TYPE >( 1, 2*n );
$TEMPLATE[end]
