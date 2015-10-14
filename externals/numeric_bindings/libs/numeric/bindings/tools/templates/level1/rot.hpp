$TEMPLATE[rot.all.remove_argument_value_type_prepend]
X,Y
$TEMPLATE[rot.all.arguments]
 N         (input) INTEGER
           The length of array X
 INCX      (input) INTEGER
           The increment of X
 INCY      (input) INTEGER
           The increment of Y
 C        (input) DATATYPE variable alpha
 S        (input) DATATYPE variable alpha
 X        (input/output) DATATYPE array of length (N)
 Y        (input/output) DATATYPE array of length (N)
$TEMPLATE[rot.all.X.io]
input,output
$TEMPLATE[rot.all.Y.io]
input,output
$TEMPLATE[rot.all.INCX.trait]
stride,X
$TEMPLATE[rot.all.INCY.trait]
stride,Y
$TEMPLATE[rot.all.N.trait]
size,X
$TEMPLATE[end]
