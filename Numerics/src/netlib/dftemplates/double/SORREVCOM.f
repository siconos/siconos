*
      SUBROUTINE SORREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                     NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      DOUBLE PRECISION   RESID
      INTEGER            NDX1, NDX2
      DOUBLE PRECISION   SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( * ), X( * ), WORK( LDW,* )
*     .. 
*
*  Purpose
*  =======
*
*  SOR solves the linear system Ax = b using the Successive 
*  Over-Relaxation iterative method. 
*  The matrix splitting is formed by copying the strict upper 
*  triangular portion of A onto matrix N, stored in WORK. Matrix M 
*  is the lower triangular portion of A.
*  On exit, matrix A and right hand side b are reset to their 
*  original form.
*
*  Relative error measured: norm( X - X_1 ) / norm( X ).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (input/workspace) DOUBLE PRECISION array, dimension (N*(N+3)).
*          The relaxation parameter, OMEGA, should be input in WORK(1).
*          The amount of workspace can be significantly reduced (to 2*N)
*          by customizing the matrix-vector product and backsolve.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION         
*          On input, the allowable convergence measure for
*          norm( x - x_1 ) / norm( x ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*            -5: Erroneous NDX1/NDX2 in INIT call.
*            -6: Erroneous RLBL.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -4: Relaxation parameter OMEGA not in interval (0,2).
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) DOUBLE PRECISION.
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS:   DAXPY, DCOPY, DNRM2
*     ==========================================================
*
*     .. External Routines ..
      EXTERNAL           DAXPY, DCOPY, DNRM2
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER        ( ONE = 1.0D+0 , ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            MM, X1, TEMP, MAXIT, NEED1, NEED2
      DOUBLE PRECISION   OMEGA, TOL, BNRM2, DNRM2
*     ..
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      X1   = 1
      TEMP = 2
      MM   = 3
*
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((X1 - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((TEMP - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((MM - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((X1 - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((TEMP - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((MM - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set relaxation parameter.
*
      OMEGA = WORK( 1,1 )
      IF ( OMEGA .EQ. ZERO ) OMEGA = ONE
*
*     Compute initial residual for ( convergence criteria ).
*
      CALL DCOPY( N, B, 1, WORK(1,X1), 1 )
      IF ( DNRM2( N, X, 1 ) .NE. ZERO ) THEN
*********CALL MATVEC( -ONE, X, ONE, WORK(1,X1) )
*
         NDX1 = -1
         NDX2 = ((X1   - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = -ONE
         SCLR2 = ONE
         RLBL = 2
         IJOB = 1
         RETURN
*
*****************
 2       CONTINUE
*****************
*
         IF ( DNRM2( N, WORK(1,X1), 1 ).LT.TOL ) GO TO 30
      ENDIF
      BNRM2 = DNRM2( N, B, 1 )
      IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
*
*     Matrix A is set to N. WORK(1:N,1:N) is set to MM.
*
      CALL MATSPLIT( OMEGA, B, WORK(1,MM), LDW,'SOR','SPLIT')
*
      ITER = 0
*
   10 CONTINUE
*
*     Perform SOR iteration
*
      ITER = ITER + 1
*
*        Save the current approximation to X in X1,
*
         CALL DCOPY( N, X, 1, WORK( 1,X1 ), 1 )
*
*        Apply iteration; result is updated approximation vector x
*
         CALL DCOPY( N, B, 1, WORK( 1,TEMP ), 1 )
*********CALL MATVEC( ONE, X, ONE, WORK( 1,TEMP ) )
*
         NDX1 = -1
         NDX2 = ((TEMP - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ONE
         RLBL = 3
         IJOB = 1
         RETURN
*
*****************
 3       CONTINUE
*****************
*
         CALL DCOPY( N, WORK( 1,TEMP ), 1, X, 1 )
*********CALL BACKSOLVE( N, WORK( 1,MM ), LDW, X )
*
         NDX1 = ((MM   - 1) * LDW) + 1
         NDX2 = -1
*        Prepare for return & return
         RLBL = 4
         IJOB = 2
         RETURN
*
*****************
 4       CONTINUE
*****************
*
*        Compute error and check for acceptable convergence.
*
         CALL DAXPY( N, -ONE, X, 1, WORK( 1,X1 ), 1 )
*
*********RESID = DNRM2( N, WORK( 1,X1 ), 1 ) / 
*****$           DNRM2( N, X, 1 )
*********IF ( RESID.LE.TOL ) GO TO 30
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 5
         IJOB = 3
         RETURN
*
*****************
 5       CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 30
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 20
         ENDIF
*
         GO TO 10
*
   20 CONTINUE
*
*     Iteration fails
*     restore B to original form
*
      CALL MATSPLIT( OMEGA, B, WORK(1,MM), LDW,'SOR','RESTORE')
*
      RLBL = -1
      IJOB = -1
*
      RETURN
*
 30   CONTINUE
*
*     Iteration successful; restore A and B to original form,
*     compute residual norm, and return
*
      CALL MATSPLIT( OMEGA, B, WORK(1,MM), LDW,'SOR','RESTORE')
*
      INFO = 0
      RLBL = -1
      IJOB = -1
*
      RETURN
*
*     End of SORREVCOM
*
      END
