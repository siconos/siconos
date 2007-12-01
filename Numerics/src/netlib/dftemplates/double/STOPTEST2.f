*
      SUBROUTINE STOPTEST2( N, R, B, BNRM2, RESID, TOL, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            N, INFO
      DOUBLE PRECISION   RESID, TOL, BNRM2
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   R( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  Computes the stopping criterion 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  INFO    (output) INTEGER
*          On exit, 1/0 depending on whether stopping criterion
*          was met or not.
*
*  BNRM2   (input/output) DOUBLE PRECISION.
*          On first time entry, will be -1.0.
*          On first time exit will contain norm2(B)
*          On all subsequent entry/exit's unchanged.
*
*  RESID   (output) DOUBLE PRECISION.
*          On exit, the computed stopping measure.
*
*  TOL     (input) DOUBLE PRECISION.
*          On input, the allowable convergence measure.
*
*  R       (input) DOUBLE PRECISION array, dimension N.
*          On entry, the residual.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  BLAS CALLS:   DNRM2
*  ============================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. External Routines ..
      DOUBLE PRECISION DNRM2
      EXTERNAL         DNRM2
*     ..
*     .. Executable Statements ..
*
      IF( INFO.EQ.-1 ) THEN
         BNRM2 = DNRM2( N, B, 1 )
         IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
      ENDIF
*
      RESID = DNRM2( N, R, 1 ) / BNRM2
*
      INFO = 0
      IF ( RESID.LE.TOL )
     $     INFO = 1
*
      RETURN
*
*     End of STOPTEST2
*
      END
