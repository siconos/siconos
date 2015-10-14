$TEMPLATE[lange.all.min_size_work.args]
NORM,M
$TEMPLATE[lange.all.min_size_work]
if ( norm == 'I' )
    return std::max< $INTEGER_TYPE >( 1, m );
else
    return 1;
$TEMPLATE[end]
