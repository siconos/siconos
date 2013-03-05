$TEMPLATE[geesx.all.min_size_bwork.args]
N, SORT
$TEMPLATE[geesx.all.min_size_bwork]
if ( sort == 'N' )
    return 0;
else
    return n;
$TEMPLATE[geesx.real.min_size_iwork.args]
N, SENSE
$TEMPLATE[geesx.real.min_size_iwork]
if ( sense == 'N' || sense == 'E' )
    return 1;
else
    return std::max< $INTEGER_TYPE >( 1, n*n/4 );
$TEMPLATE[geesx.all.min_size_work.args]
N, SENSE
$TEMPLATE[geesx.real.min_size_work]
if ( sense == 'N' )
    return std::max< $INTEGER_TYPE >( 1, 3*n );
else
    return std::max< $INTEGER_TYPE >( 1, n+n*n/2 );
$TEMPLATE[geesx.complex.min_size_work]
if ( sense == 'N' )
    return std::max< $INTEGER_TYPE >( 1, 2*n );
else
    return std::max< $INTEGER_TYPE >( 1, n*n/2 );
$TEMPLATE[end]
