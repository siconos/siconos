$TEMPLATE[trsna.all.WORK.type]
vector
$TEMPLATE[trsna.all.fixme]
LDWORK isn't handled correctly yet,
because trsna uses a matrix as workspace instead of a vector.
$TEMPLATE[trsna.all.extra_variables]
LDWORK
$TEMPLATE[trsna.all.LDWORK.init]
$INTEGER_TYPE ldwork = std::max< $INTEGER_TYPE >( 1, (job == 'V' || job == 'B') ? n : 1 );
$TEMPLATE[trsna.all.min_size_work.args]
JOB,LDWORK,N
$TEMPLATE[trsna.all.min_size_work]
return job == 'E' ? 1 : ldwork * (n+6);
$TEMPLATE[trsna.real.min_size_iwork.args]
JOB,N
$TEMPLATE[trsna.real.min_size_iwork]
return std::max< $INTEGER_TYPE >( 1, job == 'E' ? 1 : 2 * (n-1));
$TEMPLATE[trsna.complex.min_size_rwork.args]
JOB,N
$TEMPLATE[trsna.complex.min_size_rwork]
return std::max< $INTEGER_TYPE >( 1, job == 'E' ? 1 : n);
$TEMPLATE[end]
