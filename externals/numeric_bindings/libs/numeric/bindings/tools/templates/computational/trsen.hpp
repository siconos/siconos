$TEMPLATE[trsen.all.fixme]
M should be computed from select, because it is an output parameter.
$TEMPLATE[trsen.all.min_size_work.args]
JOB,N,M
$TEMPLATE[trsen.real.min_size_work]
if ( job == 'N' )
    return std::max< $INTEGER_TYPE >(1, n);
else if ( job == 'E' )
    return std::max< $INTEGER_TYPE >(1, m*(n-m));
else // if ( job == 'V' || job == 'B' )
    return std::max< $INTEGER_TYPE >(1, 2*m*(n-m));
$TEMPLATE[trsen.complex.min_size_work]
if ( job == 'N' )
    return 1;
else if ( job == 'E' )
    return std::max< $INTEGER_TYPE >(1, m*(n-m));
else // if ( job == 'V' || job == 'B' )
    return std::max< $INTEGER_TYPE >(1, 2*m*(n-m));
$TEMPLATE[trsen.real.min_size_iwork.args]
JOB,N,M
$TEMPLATE[trsen.real.min_size_iwork]
if ( job == 'N' || job == 'E' )
    return 1;
else // if ( job == 'V' || job == 'B' )
    return std::max< $INTEGER_TYPE >(1, m*(n-m));
$TEMPLATE[end]
