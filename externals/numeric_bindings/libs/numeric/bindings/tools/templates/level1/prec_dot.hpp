$TEMPLATE[prec_dot.all.remove_argument_value_type_prepend]
X,Y
$TEMPLATE[prec_dot.all.arguments]
 N         (input) INTEGER
           The length of array X
 INCX      (input) INTEGER
           The increment of X
 INCY      (input) INTEGER
           The increment of Y
 X         (input) DATATYPE array of length (N)
 Y         (input) DATATYPE array of length (N)
$TEMPLATE[prec_dot.all.cblas_alias]
SY,Y
SX,X
$TEMPLATE[prec_dot.all.N.trait]
size,X
$TEMPLATE[prec_dot.all.INCX.trait]
stride,X
$TEMPLATE[prec_dot.all.INCY.trait]
stride,Y
$TEMPLATE[prec_dot.all.level1_result_type]
double
$TEMPLATE[end]
