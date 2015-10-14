$TEMPLATE[ggevx.all.min_size_bwork.args]
SENSE, N
$TEMPLATE[ggevx.all.min_size_bwork]
if ( sense == 'N' )
  return 0;
else
  return n;
$TEMPLATE[ggevx.real.min_size_iwork.args]
SENSE, N
$TEMPLATE[ggevx.real.min_size_iwork]
if ( sense == 'E' )
  return 0;
else
  return n+6;
$TEMPLATE[ggevx.complex.min_size_iwork.args]
SENSE, N
$TEMPLATE[ggevx.complex.min_size_iwork]
if ( sense == 'E' )
  return 0;
else
  return n+2;
$TEMPLATE[ggevx.complex.min_size_rwork.args]
BALANC, N
$TEMPLATE[ggevx.complex.min_size_rwork]
if ( balanc == 'S' || balanc == 'B' )
    return std::max< $INTEGER_TYPE >( 1, 6*n );
else
    return std::max< $INTEGER_TYPE >( 1, 2*n );
$TEMPLATE[ggevx.real.min_size_work.args]
BALANC,JOBVL,JOBVR,SENSE,N
$TEMPLATE[ggevx.real.min_size_work]
if ( balanc == 'S' || balanc == 'B' || jobvl == 'V' || jobvr == 'V' )
    return std::max< $INTEGER_TYPE >( 1, 6*n );
if ( sense == 'E' )
    return std::max< $INTEGER_TYPE >( 1, 10*n );
if ( sense == 'V' || sense == 'B' )
    return 2*n*n + 8*n + 16;
return std::max< $INTEGER_TYPE >( 1, 2*n );
$TEMPLATE[ggevx.complex.min_size_work.args]
SENSE, N
$TEMPLATE[ggevx.complex.min_size_work]
if ( sense == 'N' )
    return std::max< $INTEGER_TYPE >( 1, 2*n );
else {
    if ( sense == 'E' )
        return std::max< $INTEGER_TYPE >( 1, 4*n );
    else
        return std::max< $INTEGER_TYPE >( 1, 2*n*n+2*n );
}
$TEMPLATE[end]
