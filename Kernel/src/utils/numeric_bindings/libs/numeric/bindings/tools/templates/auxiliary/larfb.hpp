$TEMPLATE[larfb.all.TRANS.trait_of]
V
$TEMPLATE[larfb.all.WORK.type]
vector
$TEMPLATE[larfb.all.fixme]
LDWORK isn't handled correctly yet,
because larfb uses a matrix as workspace instead of a vector.
$TEMPLATE[larfb.all.extra_variables]
LDWORK
$TEMPLATE[larfb.all.LDWORK.init]
$INTEGER_TYPE ldwork = std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[larfb.all.min_size_work.args]
LDWORK,K
$TEMPLATE[larfb.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[larfb.all.min_size_work]
return ldwork * k;
$TEMPLATE[end]
