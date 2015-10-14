$TEMPLATE[gejsv.all.min_size_iwork.args]
M,N
$TEMPLATE[gejsv.all.min_size_iwork]
return m+3*n;
$TEMPLATE[gejsv.all.min_size_work.args]
JOBA,JOBU,JOBV,M,N
$TEMPLATE[gejsv.all.min_size_work]
if ( jobu == 'N' && jobv == 'N' ) {
    if ( joba != 'E' && joba != 'G' )
        return std::max< $INTEGER_TYPE >( std::max< $INTEGER_TYPE >( 2*m+n, 4*n+1), 7 );
    else
        return std::max< $INTEGER_TYPE >( std::max< $INTEGER_TYPE >( 2*m+n, n*n+4*n), 7 );
} else if ( jobu == 'N' || jobu == 'W' || jobv == 'N' || jobv == 'W' ) {
        return std::max< $INTEGER_TYPE >( 2*n+m, 7);
} else {
    if ( jobv != 'J' )
        return 6*n+2*n*n;
    else
        return std::max< $INTEGER_TYPE >( m+3*n+n*n, 7);
}
$TEMPLATE[gejsv.all.A.io]
input;output
$TEMPLATE[gejsv.all.SVA.io]
output
$TEMPLATE[gejsv.all.U.io]
output
$TEMPLATE[gejsv.all.V.io]
output
$TEMPLATE[end]
