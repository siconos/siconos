$TEMPLATE[tgsna.all.min_size_work.args]
JOB,N
$TEMPLATE[tgsna.real.min_size_work]
if ( job == 'V' || job == 'B' )
    return std::max< $INTEGER_TYPE >(1, 2*n*(n+2)+16);
else
    return std::max< $INTEGER_TYPE >(1, n);
$TEMPLATE[tgsna.complex.min_size_work]
if ( job == 'V' || job == 'B' )
    return std::max< $INTEGER_TYPE >(1, 2*n*n);
else
    return std::max< $INTEGER_TYPE >(1, n);
$TEMPLATE[tgsna.all.min_size_iwork.args]
JOB,N
$TEMPLATE[tgsna.real.min_size_iwork]
if ( job == 'E')
    return 1;
else
    return n+6;
$TEMPLATE[tgsna.complex.min_size_iwork]
if ( job == 'E')
    return 1;
else
    return n+2;
$TEMPLATE[end]
