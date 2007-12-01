*
      SUBROUTINE QMR( N, B, X, WORK, LDW, ITER, RESID, MATVEC, 
     $                MATVECTRANS, PSOLVEQ, PSOLVETRANSQ, INFO )
*
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
      EXTERNAL           MATVEC, MATVECTRANS, PSOLVEQ, PSOLVETRANSQ
*
*  Purpose
*  =======
*
*  QMR Method solves the linear system Ax = b using the
*  Quasi-Minimal Residual iterative method with preconditioning.
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
*  X      (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess; on exit, the iterated solution.
*
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,11).
*          Workspace for residual, direction vector, etc.
*          Note that W and WTLD, Y and YTLD, and Z and ZTLD share
*          workspace.
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
*  MATVECTRANS  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A'*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A' is the tranpose of a matrix A. Vector x must remain
*          unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVECTRANS( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVEQ  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. As QMR uses left
*          and right preconditioning and the preconditioners are in
*          common, we must specify in the call which to use. Vector b
*          must remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVEQ( X, B, 'LEFT' )
*
*          The preconditioner is passed into the routine in a common block.
*
*  PSOLVETRANSQ  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M'*x = b,
*
*          where x and y are vectors, and M' is the tranpose of a 
*          matrix M. As QMR uses left and right preconditioning and
*          the preconditioners are in common, we must specify in the 
*          call which to use. Vector b must remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVETRANSQ( X, B, 'LEFT' )
*
*          The preconditioner is passed into the routine in a common block.
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
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO   < BREAKTOL: RHO and RTLD have become
*                                         orthogonal.
*                  -11: BETA  < BREAKTOL: EPS too small in relation to DELTA.
*                                         Convergence has stalled.
*                  -12: GAMMA < BREAKTOL: THETA too large. 
*                                         Convergence has stalled.
*                  -13: DELTA < BREAKTOL: Y and Z have become
*                                         orthogonal.
*                  -14: EPS   < BREAKTOL: Q and PTLD have become
*                                         orthogonal.
*                  -15: XI    < BREAKTOL: Z too small. Convergence has stalled.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  BLAS CALLS:   DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*  ==============================================================
*
*This variable used to communicate requests between QMR() and 
* QMRREVCOM()
*QMR -> QMRREVCOM: 1 = init, 
*                  2 = use saved state to resume flow.
*QMRREVCOM -> QMR: -1 = done, return to main, 
*                   1 = matvec using SCLR1/2, NDX1/2 
*                   2 = transpose-matvec using SCLR1/2, NDX1/2 
*                   3 = left prec solve using NDX1/2
*                   4 = right prec solve using NDX1/2
*                   5 = left prec transpose-solve using NDX1/2
*                   6 = right transpose-solve using NDX1/2
      INTEGER          IJOB
      LOGICAL          FTFLG
*     Arg/Result indices into WORK[].
      INTEGER NDX1, NDX2
*     Scalars passed from QMRREVCOM to QMR.
      DOUBLE PRECISION SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      DOUBLE PRECISION TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         QMRREVCOM, STOPTEST2
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
*     1 == R; 2 == D; 3 == P; 4 == PTLD; 5 == Q; 6 == S; 7 == V;
*     8 == VTLD; 9 == W; 10 == WTLD; 11 == Y; 12 == YTLD; 13 == Z;
*     14 == ZTLD; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

          CALL QMRREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $                   NDX1, NDX2, SCLR1, SCLR2, IJOB)

*         On a return from REVCOM we use the table (BiCGREVCOM -> BiCG)
*         to decode IJOB.
          IF (IJOB .eq. -1) THEN
*           revcom wants to terminate, so do it.
            GOTO 2
          ELSEIF (IJOB .eq. 1) THEN
*           call matvec.
            CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
          ELSEIF (IJOB .eq. 2) THEN
*           call transpose-matvec.
            CALL MATVECTRANS(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
          ELSEIF (IJOB .eq. 3) THEN
*           call left prec solve.
            CALL PSOLVEQ(WORK(NDX1), WORK(NDX2), 'LEFT')
          ELSEIF (IJOB .eq. 4) THEN
*           call right prec solve.
            CALL PSOLVEQ(WORK(NDX1), WORK(NDX2), 'RIGHT')
          ELSEIF (IJOB .eq. 5) THEN
*           call left prec transpose-solve.
            CALL PSOLVETRANSQ(WORK(NDX1), WORK(NDX2),'LEFT')
          ELSEIF (IJOB .eq. 6) THEN
*           call right prec transpose-solve.
            CALL PSOLVETRANSQ(WORK(NDX1), WORK(NDX2),'RIGHT')
          ELSEIF (IJOB .eq. 7) THEN
*           call matvec with X.
            CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
         ELSEIF (IJOB .EQ. 8) THEN
*           do stopping test 2
*           if first time, set INFO so that BNRM2 is computed.
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
*
      RETURN
*
*     End of QMR
*
      END
