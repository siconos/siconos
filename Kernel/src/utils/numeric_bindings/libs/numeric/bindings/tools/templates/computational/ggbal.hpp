$TEMPLATE[ggbal.all.min_size_work.args]
JOB,N
$TEMPLATE[ggbal.all.min_size_work]
if ( job == 'S' || job == 'B' )
    return std::max< $INTEGER_TYPE >(1, 6*n);
else // if ( job == 'N' || job == 'P' )
    return 1;
$TEMPLATE[end]
