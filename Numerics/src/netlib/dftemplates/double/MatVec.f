*
      SUBROUTINE MATVEC( ALPHA, X, BETA, Y )
*
*     This MatVec routine assumes the matrix is in dense format,
*     and uses the BLAS DGEMV.
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
*     ..
*     .. Common Blocks ..
*     MAXDIM2 = MAXDIM*MAXDIM.
*
      INTEGER             MAXDIM, MAXDIM2
      PARAMETER         ( MAXDIM = 200, MAXDIM2 = 40000 )
*
      INTEGER             N, LDA
      DOUBLE PRECISION    A, M
*
      COMMON            / SYSTEM / A( MAXDIM2 ), M( MAXDIM ),
     $                  / MATDIM / N, LDA
*     ..
      EXTERNAL           DGEMV
*
      CALL DGEMV( 'NOTRANSPOSE', N, N, ALPHA, A, LDA, X, 1, BETA, Y, 1 )
*
      RETURN
*
      END
*
*     =================================================
      SUBROUTINE MATVECTRANS( ALPHA, X, BETA, Y )
*
*     This MatVec routine assumes the matrix is in dense format,
*     and uses the BLAS DGEMV.
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
*     ..
*     .. Common Blocks ..
*     MAXDIM2 = MAXDIM*MAXDIM.
*
      INTEGER             MAXDIM, MAXDIM2
      PARAMETER         ( MAXDIM = 200, MAXDIM2 = 40000 )
*
      INTEGER             N, LDA
      DOUBLE PRECISION    A, M
*
      COMMON            / SYSTEM / A( MAXDIM2 ), M( MAXDIM ),
     $                  / MATDIM / N, LDA
*     ..
      EXTERNAL           DGEMV
*
      CALL DGEMV( 'TRANSPOSE', N, N, ALPHA, A, LDA, X, 1, BETA, Y, 1 )
*
      RETURN
*
      END
