$TEMPLATE[template_lapack_solve_pivot0]
    // high-level solve typedefs and functions
    typedef boost::mpl::bool_<false> has_pivot;

    template< typename MatrixA, typename MatrixB, typename VectorP >
    static void solve( MatrixA& A, MatrixB& B, VectorP const&, $LIBRARY_INT_TYPE& info ) {
        invoke( $KEYWORDS );
    }
$TEMPLATE[template_lapack_solve_pivot1]
    // high-level solve typedefs and functions
    typedef boost::mpl::bool_<true> has_pivot;

    template< typename MatrixA, typename MatrixB, typename VectorP >
    static void solve( MatrixA& A, MatrixB& B, VectorP& pivot, $LIBRARY_INT_TYPE& info ) {
        invoke( $KEYWORDS );
    }

    template< typename MatrixA, typename MatrixB, typename VectorP >
    static void solve( MatrixA& A, MatrixB& B, VectorP const&, $LIBRARY_INT_TYPE& info ) {
        traits::detail::array< $LIBRARY_INT_TYPE > pivot( size_column(A) );
        invoke( $KEYWORDS );
    }
$TEMPLATE[end]
