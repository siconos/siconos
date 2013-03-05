$TEMPLATE[tgsen.all.fixme]
M should be computed from select, because it is an output parameter.
$TEMPLATE[tgsen.all.min_size_work.args]
IJOB,N,M
$TEMPLATE[tgsen.real.min_size_work]
if ( ijob == 1 || ijob == 2 || ijob == 4 )
    return std::max< $INTEGER_TYPE >(4*n+16, 2*m*(n-m));
else if ( ijob == 3 || ijob == 5 )
    return std::max< $INTEGER_TYPE >(4*n+16, 4*m*(n-m));
else // ijob == 0
    return std::max< std::ptrdiff_t >(1, 4*n+16);
$TEMPLATE[tgsen.complex.min_size_work]
if ( ijob == 1 || ijob == 2 || ijob == 4 )
    return std::max< $INTEGER_TYPE >(1, 2*m*(n-m));
else if ( ijob == 3 || ijob == 5 )
    return std::max< $INTEGER_TYPE >(1, 4*m*(n-m));
else // ijob == 0
    return 1;
$TEMPLATE[tgsen.all.min_size_iwork.args]
IJOB,N,M
$TEMPLATE[tgsen.real.min_size_iwork]
if ( ijob == 1 || ijob == 2 || ijob == 4 )
    return std::max< $INTEGER_TYPE >(1, n+6);
else if ( ijob == 3 || ijob == 5 )
    return std::max< $INTEGER_TYPE >(2*m*(n-m), n+6);
else // ijob == 0
    return 1;
$TEMPLATE[tgsen.complex.min_size_iwork]
if ( ijob == 1 || ijob == 2 || ijob == 4 )
    return std::max< $INTEGER_TYPE >(1, n+2);
else if ( ijob == 3 || ijob == 5 )
    return std::max< $INTEGER_TYPE >(2*m*(n-m), n+2);
else // ijob == 0
    return 1;
$TEMPLATE[end]
