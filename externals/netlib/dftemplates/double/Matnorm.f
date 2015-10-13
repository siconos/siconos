*     ===========================================================
      DOUBLE PRECISION FUNCTION MATNORM( N, A, LDA )
*
*     Compute infinity norm of matrix A.
*
      INTEGER            N, LDA, I, J
      DOUBLE PRECISION   ROWSUM, ZERO, TEMP, A( LDA,* )
      PARAMETER        ( ZERO = 0.0D+0 )
*
      TEMP = ZERO
      DO 20 I = 1, N
         ROWSUM = ZERO
         DO 10 J = 1, N
            ROWSUM = ROWSUM + ABS( A( I,J ) )
   10    CONTINUE
         TEMP = MAX( ROWSUM, TEMP )
   20 CONTINUE
*
      MATNORM = TEMP
*
      RETURN
*
      END
