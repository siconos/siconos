$TEMPLATE[ggesx.all.min_size_iwork.args]
N, SENSE
$TEMPLATE[ggesx.real.min_size_iwork]
if ( sense == 'N' )
    return 1;
else
    return std::max< $INTEGER_TYPE >( 1, n+6 );
$TEMPLATE[ggesx.complex.min_size_iwork]
if ( sense == 'N' )
    return 1;
else
    return std::max< $INTEGER_TYPE >( 1, n+2 );
$TEMPLATE[ggesx.all.min_size_bwork.args]
N, SORT
$TEMPLATE[ggesx.all.min_size_bwork]
if ( sort == 'N' )
    return 0;
else
    return n;
$TEMPLATE[ggesx.all.min_size_work.args]
N, SENSE
$TEMPLATE[ggesx.real.min_size_work]
if ( n == 0 )
    return 1;
if ( sense == 'N' )
    return std::max< $INTEGER_TYPE >( 8*n, 6*n+16 );
else
    return std::max< $INTEGER_TYPE >( 8*n, std::max< $INTEGER_TYPE >( 6*n+16, n*n/2 ));
$TEMPLATE[ggesx.complex.min_size_work]
if ( sense == 'N' )
    return std::max< $INTEGER_TYPE >( 1, 2*n );
else
    return std::max< $INTEGER_TYPE >( 1, std::max< $INTEGER_TYPE >( 2*n, n*n/2 ) );
$TEMPLATE[end]
