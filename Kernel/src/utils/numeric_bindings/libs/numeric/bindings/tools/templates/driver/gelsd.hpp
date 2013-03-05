$TEMPLATE[gelsd.includes]
#include <boost/numeric/bindings/lapack/auxiliary/ilaenv.hpp>
$TEMPLATE[gelsd.complex.min_size_rwork.args]
MINMN,SMLSIZ,NLVL,NRHS
$TEMPLATE[gelsd.all.extra_variables]
MINMN,SMLSIZ,NLVL
$TEMPLATE[gelsd.all.extra_opt_variables]
MINMN,NLVL
$TEMPLATE[gelsd.all.MINMN.init]
$INTEGER_TYPE minmn = std::min< $INTEGER_TYPE >( size_row(a), size_column(a) );
$TEMPLATE[gelsd.all.SMLSIZ.init]
$INTEGER_TYPE smlsiz = ilaenv(9, "GELSD", "");
$TEMPLATE[gelsd.all.NLVL.init]
$INTEGER_TYPE nlvl = std::max< $INTEGER_TYPE >( static_cast<$INTEGER_TYPE>(std::log(static_cast<real_type>(minmn)/static_cast<real_type>(smlsiz+1))/std::log(2.0)) + 1, 0 );
$TEMPLATE[gelsd.complex.min_size_rwork]
$INTEGER_TYPE smlsiz_plus_one = smlsiz + 1;
return std::max< $INTEGER_TYPE >( 1, 10*minmn + 2*minmn*smlsiz + 8*minmn*nlvl + 3*smlsiz*nrhs + smlsiz_plus_one * smlsiz_plus_one );
$TEMPLATE[gelsd.complex.min_size_work.args]
N, MINMN, NRHS
$TEMPLATE[gelsd.complex.min_size_work]
return std::max< $INTEGER_TYPE >( 1, 2*minmn + std::max< $INTEGER_TYPE >( n, minmn*nrhs ) );
$TEMPLATE[gelsd.all.min_size_iwork.args]
MINMN,NLVL
$TEMPLATE[gelsd.all.min_size_iwork]
return std::max< $INTEGER_TYPE >( 1, 3*minmn*nlvl + 11*minmn );
$TEMPLATE[gelsd.real.min_size_work.args]
MINMN,SMLSIZ, NLVL, NRHS
$TEMPLATE[gelsd.real.min_size_work]
$INTEGER_TYPE smlsiz_plus_one = smlsiz + 1;
return std::max< $INTEGER_TYPE >( 1, 12*minmn + 2*minmn*smlsiz + 8*minmn*nlvl + minmn*nrhs + smlsiz_plus_one * smlsiz_plus_one );
$TEMPLATE[gelsd.real.A.io]
input;output
$TEMPLATE[end]
