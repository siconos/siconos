$TEMPLATE[bdsqr.all.min_size_work.args]
N,NCVT,NRU,NCC
$TEMPLATE[bdsqr.real.min_size_work]
if ( ncvt == 0 && nru == 0 && ncc == 0 )
    return 2*n;
else
    return std::max< $INTEGER_TYPE >(1, 4*n);
$TEMPLATE[bdsqr.complex.min_size_rwork.args]
N,NCVT,NRU,NCC
$TEMPLATE[bdsqr.complex.min_size_rwork]
if ( ncvt == 0 && nru == 0 && ncc == 0 )
    return 2*n;
else
    return std::max< $INTEGER_TYPE >(1, 4*n-4);
$TEMPLATE[end]
