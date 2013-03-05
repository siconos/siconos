$TEMPLATE[lalsd.real.min_size_work.args]
N,SMLSIZ, NLVL, NRHS
$TEMPLATE[lalsd.real.min_size_work]
$INTEGER_TYPE smlsiz_plus_one = smlsiz + 1;
return 9*n + 2*n*smlsiz + 8*n*nlvl + n*nrhs + smlsiz_plus_one * smlsiz_plus_one;
$TEMPLATE[lalsd.real.min_size_iwork.args]
N,NLVL
$TEMPLATE[lalsd.real.min_size_iwork]
return 3*n*nlvl + 11*n;
$TEMPLATE[lalsd.complex.min_size_rwork.args]
N,SMLSIZ, NLVL, NRHS
$TEMPLATE[lalsd.complex.min_size_rwork]
$INTEGER_TYPE smlsiz_plus_one = smlsiz + 1;
return 9*n + 2*n*smlsiz + 8*n*nlvl + 3*smlsiz*nrhs + smlsiz_plus_one * smlsiz_plus_one;
$TEMPLATE[lalsd.all.extra_variables]
NLVL
$TEMPLATE[lalsd.complex.NLVL.init]
$INTEGER_TYPE nlvl = std::max< $INTEGER_TYPE >( 0, static_cast<$INTEGER_TYPE>(
    std::log(static_cast<real_type>(n)/
    static_cast<real_type>(smlsiz+1)) /
    std::log(static_cast<real_type>(2.))) + 1 );
$TEMPLATE[lalsd.real.NLVL.init]
$INTEGER_TYPE nlvl = std::max< $INTEGER_TYPE >( 0, static_cast<$INTEGER_TYPE>(
    std::log(static_cast<real_type>(n)/static_cast<real_type>(smlsiz+1)) /
    std::log(static_cast<real_type>(2.)) ) + 1 );
$TEMPLATE[end]
