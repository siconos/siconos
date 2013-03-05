$TEMPLATE[langb.all.min_size_work.args]
NORM,N
$TEMPLATE[langb.all.min_size_work]
if ( norm == 'I' )
    return std::max< $INTEGER_TYPE >( 1, n );
else
    return 1;
$TEMPLATE[end]
