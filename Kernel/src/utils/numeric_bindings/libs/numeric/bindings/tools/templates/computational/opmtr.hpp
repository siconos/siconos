$TEMPLATE[opmtr.all.TRANS.trait_of]
AP
$TEMPLATE[opmtr.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[opmtr.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[opmtr.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[end]
