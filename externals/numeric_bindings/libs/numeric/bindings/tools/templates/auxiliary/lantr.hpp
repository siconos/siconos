$TEMPLATE[lantr.all.UPLO.trait_of]
A
$TEMPLATE[lantr.all.min_size_work.args]
NORM,M
$TEMPLATE[lantr.all.min_size_work]
if ( norm == 'I' )
    return std::max< $INTEGER_TYPE >( 1, m );
else
    return 1;
$TEMPLATE[end]
