*
      SUBROUTINE CGS( N, B, X, WORK, LDW, ITER, RESID, MATVEC, PSOLVE,
     $                INFO )
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
      DOUBLE PRECISION   X( * ), B( * ), WORK( * )
*     ..
*     .. Function Arguments ..
      EXTERNAL           MATVEC, PSOLVE
*     ..
*
*  Purpose
*  =======
*
*  CGS solves the linear system Ax = b using the
*  Conjugate Gradient Squared iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector. The user should be warned that for
*          this particular algorithm, an initial guess close to
*          the actual solution can result in divergence.
*          On exit, the iterated solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors PHAT and QHAT, and UHAT and VHAT share 
*          the same workspace.
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
*          norm( b - A*x ) / norm( b ).
*          On ouput, the final value of this measure.
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
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit.
*          >  0: Convergence not achieved. This will be set
*                to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occured
*                during iteration.
*
*                Illegal parameter: 
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*                BREAKDOWN: If RHO become smaller than some tolerance,
*                   the program will terminate. Here we check 
*                   against tolerance BREAKTOL.
*
*                   -10: RHO < BREAKTOL: RHO and RTLD have become
*                                        orthogonal.
*
*  =============================================================
*     .. Local Scalars ..
*
*This variable used to communicate requests between CGS() and CGSREVCOM()
*CGS -> CGSREVCOM: 1 = init, 
*                2 = use saved state to resume flow.
*CGSREVCOM -> CGS: -1 = done, return to main, 
*                 1 = matvec using SCLR1/2, NDX1/2 
*                 2 = solve using NDX1/2
      INTEGER          IJOB
      LOGICAL          FTFLG
*
*     Arg/Result indices into WORK[].
      INTEGER NDX1, NDX2
*     Scalars passed from CGSREVCOM to CGS.
      DOUBLE PRECISION SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      DOUBLE PRECISION TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         CGSREVCOM, STOPTEST2
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
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == R; 2 == RTLD; 3 == P; 4 == PHAT; 5 == Q; 6 == QHAT;
*     7 == U; 8 == UHAT; 9 == VHAT; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

          CALL CGSREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $                   NDX1, NDX2, SCLR1, SCLR2, IJOB)

*         On a return from REVCOM() we use the table
*         to figure out what is reqd.
          IF (IJOB .eq. -1) THEN
*           revcom wants to terminate, so do it.
            GOTO 2
          ELSEIF (IJOB .eq. 1) THEN
*           call matvec.
            CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
          ELSEIF (IJOB .eq. 2) THEN
*           call solve.
            CALL PSOLVE(WORK(NDX1), WORK(NDX2))
          ELSEIF (IJOB .eq. 3) THEN
*           call matvec with X.
            CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
         ELSEIF (IJOB .EQ. 4) THEN
*           do stopping test 2
*           if FirstTime, set info  so that BNRM2 is computed.
            IF( FTFLG ) INFO = -1
            CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO)
            FTFLG = .FALSE.
         ENDIF
*
*        done what revcom asked, set IJOB & go back to it.
         IJOB = 2
         GOTO 1
*
*     come here to terminate
 2    CONTINUE
*
      RETURN
*
*     End of CGS
*
      END






