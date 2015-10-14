$TEMPLATE[larfx.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[larfx.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[larfx.all.min_size_work]
$INTEGER_TYPE order = bindings::detail::if_left( side, n, m );
if ( order < 11)
    return 1;
else
    return std::max< $INTEGER_TYPE >( 1, order );
$TEMPLATE[larfx.all.LDC.assert_ge]
max(1,M)
$TEMPLATE[end]
