$TEMPLATE[geevx.real.min_size_iwork.args]
SENSE,N
$TEMPLATE[geevx.real.min_size_iwork]
if ( sense == 'N' || sense == 'E' )
    return 0;
else
    return 2*n-2;
$TEMPLATE[geevx.real.min_size_work.args]
SENSE,JOBVL,JOBVR,N
$TEMPLATE[geevx.real.min_size_work]
if ( sense == 'N' || sense == 'E' ) {
    if ( jobvl =='V' || jobvr == 'V' )
        return std::max< $INTEGER_TYPE >( 1, 3*n );
    else
        return std::max< $INTEGER_TYPE >( 1, 2*n );
} else
    return std::max< $INTEGER_TYPE >( 1, n*(n+6) );
$TEMPLATE[geevx.complex.min_size_work.args]
SENSE,N
$TEMPLATE[geevx.complex.min_size_work]
if ( sense == 'N' || sense == 'E' )
    return std::max< $INTEGER_TYPE >( 1, 2*n );
else
    return std::max< $INTEGER_TYPE >( 1, n*n + 2*n );
$TEMPLATE[end]
