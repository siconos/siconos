$TEMPLATE[ormbr.all.fixme]
Declaring NQ as extra variables leads to a key error in the generator,
so the init below is pretty useless right now. It only shows how it could work in theory.
$TEMPLATE[ormbr.all.NQ.init]
$INTEGER_TYPE nq = bindings::detail::if_left( side, bindings::size_row(c), bindings::size_column(c) );
$TEMPLATE[ormbr.all.K.trait]
size,TAU
$TEMPLATE[ormbr.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[ormbr.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[ormbr.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[ormlq.all.K.trait]
size,TAU
$TEMPLATE[ormlq.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[ormlq.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[ormlq.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[ormqr.all.K.trait]
size,TAU
$TEMPLATE[ormqr.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[ormqr.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[ormqr.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[ormrq.all.K.trait]
size,TAU
$TEMPLATE[ormrq.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[ormrq.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[ormrq.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[ormrz.all.K.trait]
size,TAU
$TEMPLATE[ormrz.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[ormrz.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[ormrz.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[ormtr.all.K.trait]
size,TAU
$TEMPLATE[ormtr.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[ormtr.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[ormtr.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[ormql.all.A.io]
input;output
$TEMPLATE[ormql.all.K.trait]
size,TAU
$TEMPLATE[ormql.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[ormql.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[ormql.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[ormhr.all.K.trait]
size,TAU
$TEMPLATE[ormhr.all.min_size_work.args]
SIDE,M,N
$TEMPLATE[ormhr.includes]
#include <boost/numeric/bindings/detail/if_left.hpp>
$TEMPLATE[ormhr.all.min_size_work]
return std::max< $INTEGER_TYPE >( 1, bindings::detail::if_left( side, n, m ) );
$TEMPLATE[end]
