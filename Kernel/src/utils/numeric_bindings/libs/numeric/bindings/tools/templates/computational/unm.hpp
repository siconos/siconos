$TEMPLATE[unmbr.all.fixme]
Declaring NQ as extra variables leads to a key error in the generator,
so the init below is pretty useless right now. It only shows how it could work in theory.
$TEMPLATE[unmbr.all.NQ.init]
$INTEGER_TYPE nq = bindings::detail::if_left( side, bindings::size_row(c), bindings::size_column(c) );
$TEMPLATE[unmbr.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[unmbr.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[unmbr.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[unmhr.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[unmhr.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[unmhr.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[unmlq.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[unmlq.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[unmlq.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[unmqr.all.K.trait]
size,TAU
$TEMPLATE[unmqr.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[unmqr.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[unmqr.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[unmrq.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[unmrq.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[unmrq.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[unmrz.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[unmrz.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[unmrz.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[unmql.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[unmql.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[unmql.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[unmtr.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[unmtr.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[unmtr.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[end]
