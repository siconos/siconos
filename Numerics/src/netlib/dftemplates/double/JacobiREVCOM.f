*
      SUBROUTINE JACOBIREVCOM(N, B, X, WORK, LDW, ITER, RESID,
     $                     INFO, NDX1, NDX2, SCLR1, SCLR2, IJOB)
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
*
*  Purpose
*  =======
*
*  JACOBI solves the linear system Ax = b using the Jacobi iterative 
*  method. The matrix splitting should be accomplished before calling
*  this routine. The diagonal elements of the matrix must be passed into
*  this routine in the first column of matrix WORK. 
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
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,3)
*          Workspace for residual, direction vector, etc.
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
*  ==========================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, MAXIT, X1, MM, TEMP, NEED1, NEED2
      DOUBLE PRECISION   TOL, DNRM2
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*     ..
*     .. External Routines ..
*
      EXTERNAL           DAXPY, DCOPY, DNRM2, MATSPLIT
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
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
      MM   = 1
      X1   = 2
      TEMP = 3
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((MM - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((X1 - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((TEMP - 1) * LDW) + 1
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
            NEED2 = ((MM - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((X1 - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((TEMP - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
      ITER = 0
*
*     Form matrix splitting inv(M) and N.
*
      CALL MATSPLIT( ONE, B, WORK( 1,MM ), LDW,'JACOBI','SPLIT' )
*
   10 CONTINUE
*
*        Perform Jacobi iteration
*
         ITER = ITER + 1
*
*        Save the current approximation to X in X1.
*
         CALL DCOPY( N, X, 1, WORK( 1,X1 ),  1 )
*
*        Apply iteration; result is updated approximation vector x.
*
         CALL DCOPY( N, B, 1, WORK( 1,TEMP ), 1 )
*********CALL MATVEC( ONE, X, ONE, WORK( 1,TEMP ) )
*
         NDX1 = -1
         NDX2 = ((TEMP - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ONE
         RLBL = 2
         IJOB = 1
         RETURN
*
*****************
 2       CONTINUE
*****************
*
         DO 15 I = 1, N
            X( I ) = WORK( I,MM ) *  WORK( I,TEMP )
   15    CONTINUE
*
*        Compute error and check for acceptable convergence.
*
         CALL DAXPY( N, -ONE, X, 1, WORK( 1,X1 ), 1 )
*
*********RESID = DNRM2( N, WORK( 1,X1 ), 1 ) / DNRM2( N, X, 1 )
*********IF ( RESID.LE.TOL  ) GO TO 30
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 3
         IJOB = 2
         RETURN
*
*****************
 3       CONTINUE
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
*     Iteration fails.
*     Reconstruct matrix A.
*
      CALL MATSPLIT( ONE, B, WORK( 1,MM ), LDW, 'JACOBI','RECONSTRUCT' )
*
      RLBL = -1
      IJOB = -1
*
      RETURN
*
   30 CONTINUE 
*
*     Iteration successful. Reconstruct matrix A.
*
      CALL MATSPLIT( ONE, B, WORK( 1,MM ), LDW, 'JACOBI','RECONSTRUCT' )
*
      INFO = 0
      RLBL = -1
      IJOB = -1
*
      RETURN
*
*     End of JACOBIREVCOM
*
      END


