*
      SUBROUTINE BICG( N, B, X, WORK, LDW, ITER, RESID, MATVEC, 
     $                 MATVECTRANS, PSOLVE, PSOLVETRANS, INFO )
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
*
      EXTERNAL           MATVEC, MATVECTRANS, PSOLVE, PSOLVETRANS
*     ..
*
*  Purpose
*  =======
*
*  BiCG solves the linear system Ax = b using the
*  BiConjugate Gradient iterative method with preconditioning.
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
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,6).
*          Workspace for residual, direction vector, etc.
*          Note that Z and Q, and ZTLD and QTLD share workspace.
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
*  PSOLVETRANS  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M'*x = b,
*
*          where x and y are vectors, and M' is the tranpose of a 
*          matrix M. Vector b must remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVETRANS( X, B )
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
*                  -10: RHO < BREAKTOL: RHO and RTLD have become
*                                       orthogonal.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  ==============================================================
*
*
*This variable used to communicate requests between BiCG() and 
* BiCG_CGREVCOM()
*BiCG -> BiCGREVCOM: 1 = init, 
*                    2 = use saved state to resume flow.
*BiCGREVCOM -> BiCG: -1 = done, return to main, 
*                     1 = matvec using SCLR1/2, NDX1/2 
*                     2 = transpose-matvec using SCLR1/2, NDX1/2 
*                     3 = solve using NDX1/2
*                     4 = transpose-solve using NDX1/2
      INTEGER          IJOB
      LOGICAL          FTFLG
*
*     Arg/Result indices into WORK[].
      INTEGER          NDX1, NDX2
*     Scalars passed from BiCGREVCOM to BiCG.
      DOUBLE PRECISION SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      DOUBLE PRECISION TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         BICGREVCOM, STOPTEST2
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
*     1 == R; 2 == RTLD; 3 == Z; 4 == ZTLD; 5 == P; 6 == PTLD; 
*     7 == Q; 8 == QTLD; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

      CALL BICGREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $                NDX1, NDX2, SCLR1, SCLR2, IJOB)

*     On a return from REVCOM we use the table (BiCGREVCOM -> BiCG)
*     to decode IJOB.
      IF (IJOB .eq. -1) THEN
*        revcom wants to terminate, so do it.
         GOTO 2
      ELSEIF (IJOB .eq. 1) THEN
*        call matvec.
         CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
      ELSEIF (IJOB .eq. 2) THEN
*        call transpose-matvec.
         CALL MATVECTRANS(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
      ELSEIF (IJOB .eq. 3) THEN
*        call solve.
         CALL PSOLVE(WORK(NDX1), WORK(NDX2))
      ELSEIF (IJOB .eq. 4) THEN
*        call transpose-solve.
         CALL PSOLVETRANS(WORK(NDX1), WORK(NDX2))
      ELSEIF (IJOB .eq. 5) THEN
*        call matvec with X.
         CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
      ELSEIF (IJOB .EQ. 6) THEN
*        do stopping test 2
*        if FirstTime, set info  so that BNRM2 is computed.
         IF( FTFLG ) INFO = -1
         CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO)
         FTFLG = .FALSE.
      ENDIF
*
*     Done what revcom asked, set IJOB & go back to it.
*
      IJOB = 2
      GOTO 1
*
*     come here to terminate
 2    CONTINUE
*
      RETURN
*
*     End of BICG
*
      END
