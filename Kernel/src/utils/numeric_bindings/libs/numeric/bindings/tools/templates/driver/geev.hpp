$TEMPLATE[geev.real.min_size_work.args]
JOBVL,JOBVR,N
$TEMPLATE[geev.real.min_size_work]
if ( jobvl == 'V' || jobvr == 'V' )
    return std::max< $INTEGER_TYPE >( 1, 4*n );
else
    return std::max< $INTEGER_TYPE >( 1, 3*n );
$TEMPLATE[geev.complex.min_size_work.args]
N
$TEMPLATE[geev.complex.min_size_work]
return std::max< $INTEGER_TYPE >( 1, 2*n );
$TEMPLATE[end]
