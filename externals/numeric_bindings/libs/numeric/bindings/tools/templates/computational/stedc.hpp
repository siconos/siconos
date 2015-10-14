$TEMPLATE[stedc.all.N.trait]
size,D
$TEMPLATE[stedc.all.min_size_iwork.args]
COMPZ,N
$TEMPLATE[stedc.all.min_size_iwork]
if ( compz == 'N' || n <= 1 ) {
    return 1;
} else if ( compz == 'V' ) {
   return 6 + 6*n + 5*n*static_cast<$INTEGER_TYPE>(std::ceil(std::log(n)/std::log(2)));
} else { // compz == 'I'
   return 3 + 5*n;
}
$TEMPLATE[stedc.all.min_size_work.args]
COMPZ,N
$TEMPLATE[stedc.real.min_size_work]
if ( compz == 'N' || n <= 1 ) {
    return 1;
} else if ( compz == 'V' ) {
   return 1 + 3*n + 2*n*static_cast<$INTEGER_TYPE>(std::ceil(std::log(n)/std::log(2))) + 3*n*n;
} else { // compz == 'I'
   return 1 + 4*n + n*n;
}
$TEMPLATE[stedc.complex.min_size_work]
if ( compz == 'N' || compz == 'I' || n <= 1 ) {
    return 1;
} else { // compz == 'V'
   return n*n;
}
$TEMPLATE[stedc.complex.min_size_rwork.args]
COMPZ,N
$TEMPLATE[stedc.complex.min_size_rwork]
if ( compz == 'N' || n <= 1 ) {
    return 1;
} else if ( compz == 'V' ) {
   return 1 + 3*n + 2*n*static_cast<$INTEGER_TYPE>(std::ceil(std::log(n)/std::log(2))) + 3*n*n;
} else { // compz == 'I'
   return 1 + 4*n + 2*n*n;
}
$TEMPLATE[end]
