$TEMPLATE[dot.all.remove_argument_value_type_prepend]
X,Y
$TEMPLATE[dot.all.arguments]
 N         (input) INTEGER
           The length of array X
 INCX      (input) INTEGER
           The increment of X
 INCY      (input) INTEGER
           The increment of Y
 X         (input) DATATYPE array of length (N)
 Y         (input) DATATYPE array of length (N)
$TEMPLATE[dot.level0.gsub]
return cblas_cdotu_sub( n, x, incx, y, incy );->std::complex<float> result;
    cblas_cdotu_sub( n, x, incx, y, incy, &result );
    return result;--return cblas_zdotu_sub( n, x, incx, y, incy );->std::complex<double> result;
    cblas_zdotu_sub( n, x, incx, y, incy, &result );
    return result;--
$TEMPLATE[end]
