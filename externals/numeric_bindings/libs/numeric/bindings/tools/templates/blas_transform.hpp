$TEMPLATE[template_blas_transform]
    // high-level transform typedefs and functions
    template< typename MatrixA, typename VectorX, typename VectorY >
    static result_type transform( MatrixA& A, VectorX& x, VectorY& y, const value_type alpha, const value_type beta ) {
        invoke( $KEYWORDS );
    }
$TEMPLATE[end]
