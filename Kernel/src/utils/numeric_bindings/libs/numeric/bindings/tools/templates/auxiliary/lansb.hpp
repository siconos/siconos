$TEMPLATE[lansb.all.min_size_work.args]
NORM,N
$TEMPLATE[lansb.all.min_size_work]
if ( norm == 'I' || norm == '1' || norm == 'O' )
    return std::max< $INTEGER_TYPE >( 1, n );
else
    return 1;
$TEMPLATE[end]
