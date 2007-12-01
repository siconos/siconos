*
      SUBROUTINE SOR( N, B, X, WORK, LDW, ITER, RESID, MATVEC, 
     $                BACKSOLVE, INFO )
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
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( * ), X( * ), WORK( * )
*     .. 
*     .. Function Arguments ..
      EXTERNAL           MATVEC, BACKSOLVE
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
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  BACKSOLVE (external subroutine)
*          The user must provide a subroutine to perform the
*          linear system solve
*
*               x := M*x,
*
*          where x is a vector and M is a lower triangular matrix. 
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL BACKSOLVE( N, M, LDM, X )
*
*          The matrix is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
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
*     ==========================================================
*
*
*This variable used to communicate requests between SOR() and SORREVCOM()
*SOR -> SORREVCOM: 1 = init, 
*                2 = use saved state to resume flow.
*SORREVCOM -> SOR: -1 = done, return to main, 
*                 1 = matvec using SCLR1/2, NDX1/2 
*                 2 = solve using NDX1/2
      INTEGER          IJOB
*
*     Arg/Result indices into WORK[].
      INTEGER NDX1, NDX2
*     Scalars passed from SORREVCOM to SOR.
      DOUBLE PRECISION SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      DOUBLE PRECISION TOL, XNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         SORREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ELSE IF ( ( WORK( 1 ).LE.0 ).OR.( WORK( 1 ).GE.2 ) ) THEN
         INFO = -4
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == X1; 2 == TEMP; 3 == MM; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

          CALL SORREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $                   NDX1, NDX2, SCLR1, SCLR2, IJOB)

*         On a return from REVCOM() we use the table
*         to decode IJOB.
          IF (IJOB .eq. -1) THEN
*           revcom wants to terminate, so do it.
            GOTO 2
          ELSEIF (IJOB .eq. 1) THEN
*           call matvec with X.
            CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
          ELSEIF (IJOB .eq. 2) THEN
*           call backsolve with X
            CALL BACKSOLVE(N, WORK(NDX1), LDW, X)
         ELSEIF (IJOB .EQ. 3) THEN
*           do stopping test 2
*           XNRM2 shld be recomputed everytime.
            INFO = -1
            CALL STOPTEST2(N, WORK(NDX1), X, XNRM2, RESID, TOL, INFO)
         ENDIF
*
*     done what revcom asked, set IJOB & go back to it.
      IJOB = 2
      GOTO 1
*
*     come here to terminate
 2    CONTINUE
*
*
      RETURN
*
*     End of SOR
*
      END
