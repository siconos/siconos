$TEMPLATE[axpy.all.remove_argument_value_type_prepend]
A,X,Y
$TEMPLATE[axpy.all.arguments]
 N         (input) INTEGER
           The length of array X
 INCX      (input) INTEGER
           The increment of X
 INCY      (input) INTEGER
           The increment of Y
 A        (input) DATATYPE variable alpha
 X        (input) DATATYPE array of length (N)
 Y        (output) DATATYPE array of length (N)
$TEMPLATE[axpy.all.cblas_alias]
A,ALPHA
$TEMPLATE[axpy.friendly_name]
a times x plus y
$TEMPLATE[end]
