$TEMPLATE[dotc.all.remove_argument_value_type_prepend]
X,Y
$TEMPLATE[dotc.all.arguments]
 N         (input) INTEGER
           The length of array X
 INCX      (input) INTEGER
           The increment of X
 INCY      (input) INTEGER
           The increment of Y
 X         (input) DATATYPE
 Y         (input) DATATYPE
$TEMPLATE[dotc.friendly_name]
TODO
$TEMPLATE[dotc.level0.gsub]
return cblas_cdotc_sub( n, x, incx, y, incy );->std::complex<float> result;
    cblas_cdotc_sub( n, x, incx, y, incy, &result );
    return result;--return cblas_zdotc_sub( n, x, incx, y, incy );->std::complex<double> result;
    cblas_zdotc_sub( n, x, incx, y, incy, &result );
    return result;--
$TEMPLATE[end]
