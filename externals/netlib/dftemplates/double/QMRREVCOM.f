*
      SUBROUTINE QMRREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                     NDX1, NDX2, SCLR1, SCLR2, IJOB)
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
      INTEGER            NDX1, NDX2
      DOUBLE PRECISION   SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), B( * ), WORK( LDW,* )
*     ..
*  Purpose
*  =======
*
*  QMR Method solves the linear system Ax = b using the
*  Quasi-Minimal Residual iterative method with preconditioning.
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
*  BLAS CALLS:   DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*  ==============================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE, ZERO
      PARAMETER      ( ONE = 1.0D+0 , ZERO = 0.0D+0)
*
*     .. Local Scalars ..
      INTEGER          R, D, P, PTLD, Q, S, V, VTLD, W, WTLD, Y, YTLD,
     $                 Z, ZTLD, MAXIT, NEED1, NEED2
      DOUBLE PRECISION TOL, BETA, GAMMA, GAMMA1, DELTA, EPS, ETA, XI, 
     $                 RHO, RHO1, THETA, THETA1, BNRM2,
     $                 RHOTOL, BETATOL, GAMMATOL, DELTATOL, 
     $                 EPSTOL, XITOL, GETBREAK,
     $                 DDOT, DNRM2, C1
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*
*     ..
*     .. External Routines ..
      EXTERNAL         DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        ABS, SQRT
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
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
         IF (RLBL .eq. 8) GOTO 8
         IF (RLBL .eq. 9) GOTO 9
         IF (RLBL .eq. 10) GOTO 10
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
      R     = 1
      D     = 2
      P     = 3
      PTLD  = 4
      Q     = 5
      S     = 6
      V     = 7
      VTLD  = 8
      W     = 9
      WTLD  = 9
      Y     = 10
      YTLD  = 10
      Z     = 11
      ZTLD  = 11
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((D - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((PTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((Q - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((S - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.7 ) THEN
            NEED1 = ((V - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.8 ) THEN
            NEED1 = ((VTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.9 ) THEN
            NEED1 = ((W - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.10 ) THEN
            NEED1 = ((WTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.11 ) THEN
            NEED1 = ((Y - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.12 ) THEN
            NEED1 = ((YTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.13 ) THEN
            NEED1 = ((Z - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.14 ) THEN
            NEED1 = ((ZTLD - 1) * LDW) + 1
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
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((D - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((PTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((Q - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((S - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.7 ) THEN
            NEED2 = ((V - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.8 ) THEN
            NEED2 = ((VTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.9 ) THEN
            NEED2 = ((W - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.10 ) THEN
            NEED2 = ((WTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.11 ) THEN
            NEED2 = ((Y - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.12 ) THEN
            NEED2 = ((YTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.13 ) THEN
            NEED2 = ((Z - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.14 ) THEN
            NEED2 = ((ZTLD - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set breakdown tolerances.
*
      RHOTOL   = GETBREAK()
      BETATOL  = GETBREAK()
      GAMMATOL = GETBREAK()
      DELTATOL = GETBREAK()
      EPSTOL   = GETBREAK()
      XITOL    = GETBREAK()
*
*     Set initial residual.
*
      CALL DCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( DNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL MATVEC( -ONE, X, ZERO, WORK(1,R) )
*        Note: using D as temp
*********CALL DCOPY( N, X, 1, WORK(1,D), 1 )
         SCLR1 = -ONE
         SCLR2 = ZERO
         NDX1 = ((D - 1) * LDW) + 1
         NDX2 = ((R - 1) * LDW) + 1
         RLBL = 2
         IJOB = 7
         RETURN
*****************
 2       CONTINUE
*****************
*
         IF ( DNRM2( N, WORK(1,R), 1 ) .LT. TOL ) GO TO 30
      ENDIF
*
      BNRM2 = DNRM2( N, B, 1 )
      IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
*
      CALL DCOPY( N, WORK(1,R), 1, WORK(1,VTLD), 1 )
******CALL PSOLVEQ( WORK(1,Y), WORK(1,VTLD), 'LEFT' )
*
         NDX1 = ((Y    - 1) * LDW) + 1
         NDX2 = ((VTLD - 1) * LDW) + 1
         RLBL = 3
         IJOB = 3
         RETURN
*****************
 3       CONTINUE
*****************
*
      RHO   = DNRM2( N, WORK(1,Y), 1 )
*
      CALL DCOPY( N, WORK(1,R), 1, WORK(1,WTLD), 1 )
******CALL PSOLVETRANSQ( WORK(1,Z), WORK(1,WTLD), 'RIGHT' )
*
         NDX1 = ((Z    - 1) * LDW) + 1
         NDX2 = ((WTLD - 1) * LDW) + 1
         RLBL = 4
         IJOB = 6
         RETURN
*****************
 4       CONTINUE
*****************
*
      XI   = DNRM2( N, WORK(1,Z), 1 )
*
      GAMMA = ONE
      ETA   = -ONE
      THETA = ZERO
*
      ITER = 0
*
   40 CONTINUE
*
*     Perform Preconditioned QMR iteration.
*
         ITER = ITER + 1
*
         IF ( ( ABS( RHO ).LT.RHOTOL ).OR.( ABS( XI ).LT.XITOL ) )
     $      GO TO 25
*
         CALL DCOPY( N, WORK(1,VTLD), 1, WORK(1,V), 1 )
         CALL DSCAL( N, ONE / RHO, WORK(1,V), 1 )
         CALL DSCAL( N, ONE / RHO, WORK(1,Y), 1 )
*
         CALL DCOPY( N, WORK(1,WTLD), 1, WORK(1,W), 1 )
         CALL DSCAL( N, ONE / XI, WORK(1,W), 1 )
         CALL DSCAL( N, ONE / XI, WORK(1,Z), 1 )
*
         DELTA = DDOT( N, WORK(1,Z), 1, WORK(1,Y), 1 )
         IF ( ABS( DELTA ).LT.DELTATOL ) GO TO 25
*
*********CALL PSOLVEQ( WORK(1,YTLD), WORK(1,Y), 'RIGHT' )
*
         NDX1 = ((YTLD - 1) * LDW) + 1
         NDX2 = ((Y    - 1) * LDW) + 1
         RLBL = 5
         IJOB = 4
         RETURN
*****************
 5       CONTINUE
*****************
*
*********CALL PSOLVETRANSQ( WORK(1,ZTLD), WORK(1,Z), 'LEFT' )
*
         NDX1 = ((ZTLD - 1) * LDW) + 1
         NDX2 = ((Z    - 1) * LDW) + 1
         RLBL = 6
         IJOB = 5
         RETURN
*****************
 6       CONTINUE
*****************
*
*
         IF ( ITER.GT.1 ) THEN
            C1 = -( XI * DELTA / EPS )
            CALL DAXPY( N, C1, WORK(1,P), 1, WORK(1,YTLD), 1 )
            CALL DCOPY( N, WORK(1,YTLD), 1, WORK(1,P), 1 )
            CALL DAXPY( N, -( RHO * DELTA / EPS ), WORK(1,Q), 1, 
     $                  WORK(1,ZTLD), 1 )
            CALL DCOPY( N, WORK(1,ZTLD), 1, WORK(1,Q), 1 )
         ELSE
            CALL DCOPY( N, WORK(1,YTLD), 1, WORK(1,P), 1 )
            CALL DCOPY( N, WORK(1,ZTLD), 1, WORK(1,Q), 1 )
         ENDIF
*
*********CALL MATVEC( ONE, WORK(1,P), ZERO, WORK(1,PTLD) )
*
         SCLR1 = ONE
         SCLR2 = ZERO
         NDX1 = ((P    - 1) * LDW) + 1
         NDX2 = ((PTLD - 1) * LDW) + 1
         RLBL = 7
         IJOB = 1
         RETURN
*****************
 7       CONTINUE
*****************
*
*
         EPS = DDOT( N, WORK(1,Q), 1, WORK(1,PTLD), 1 )
         IF ( ABS( EPS ).LT.EPSTOL ) GO TO 25
*
         BETA = EPS / DELTA
         IF ( ABS( BETA ).LT.BETATOL ) GO TO 25
*
         CALL DCOPY( N, WORK(1,PTLD), 1, WORK(1,VTLD), 1 )
         CALL DAXPY( N, -BETA, WORK(1,V), 1, WORK(1,VTLD), 1 )
         CALL PSOLVEQ( WORK(1,Y), WORK(1,VTLD), 'LEFT' )
*
         RHO1 = RHO
         RHO  = DNRM2( N, WORK(1,Y), 1 )
*
         CALL DCOPY( N, WORK(1,W), 1, WORK(1,WTLD), 1 )
*********CALL MATVECTRANS( ONE, WORK(1,Q), -BETA, WORK(1,WTLD) )
*
         SCLR1 = ONE
         SCLR2 = -BETA
         NDX1 = ((Q    - 1) * LDW) + 1
         NDX2 = ((WTLD - 1) * LDW) + 1
         RLBL = 8
         IJOB = 2
         RETURN
*****************
 8       CONTINUE
*****************
*
*********CALL PSOLVETRANSQ( WORK(1,Z), WORK(1,WTLD), 'RIGHT' )
*
         NDX1 = ((Z    - 1) * LDW) + 1
         NDX2 = ((WTLD - 1) * LDW) + 1
         RLBL = 9
         IJOB = 6
         RETURN
*****************
 9       CONTINUE
*****************
*
*
         XI = DNRM2( N, WORK(1,Z), 1 )
*
         GAMMA1 = GAMMA
         THETA1 = THETA
*
         THETA  = RHO / ( GAMMA1 * ABS( BETA ) )
         GAMMA = ONE / SQRT( ONE + THETA**2 )
         IF ( ABS( GAMMA ).LT.GAMMATOL ) GO TO 25
*
         ETA = -ETA * RHO1 * GAMMA**2 / ( BETA * GAMMA1**2 )
*
         IF ( ITER.GT.1 ) THEN
            CALL DSCAL( N, ( THETA1*GAMMA )**2, WORK(1,D), 1 )
            CALL DAXPY( N, ETA, WORK(1,P), 1, WORK(1,D), 1 )
            CALL DSCAL( N, ( THETA1 * GAMMA )**2, WORK(1,S), 1 )
            CALL DAXPY( N, ETA, WORK(1,PTLD), 1, WORK(1,S), 1 )
         ELSE
            CALL DCOPY( N, WORK(1,P), 1, WORK(1,D), 1 )
            CALL DSCAL( N, ETA, WORK(1,D), 1 )
            CALL DCOPY( N, WORK(1,PTLD), 1, WORK(1,S), 1 )
            CALL DSCAL( N, ETA, WORK(1,S), 1 )
         ENDIF
*
*        Compute current solution vector x.
*
         CALL DAXPY( N, ONE, WORK(1,D), 1, X, 1 )
*
*        Compute residual vector rk, find norm,
*        then check for tolerance.
*
         CALL DAXPY( N, -ONE, WORK(1,S), 1, WORK(1,R), 1 )
*
*********RESID = DNRM2( N, WORK(1,R), 1 ) / BNRM2
*********IF ( RESID .LE. TOL  ) GO TO 30
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 10
         IJOB = 8
         RETURN
*
*****************
 10      CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 30
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 20
         ENDIF
*
         GO TO 40
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
*
      RETURN
*
   25 CONTINUE
*
*     Method breakdown.
*
      IF ( ABS( RHO ).LT.RHOTOL ) THEN
         INFO = -10
      ELSE IF ( ABS( BETA ).LT.BETATOL ) THEN
         INFO = -11
      ELSE IF ( ABS( GAMMA ).LT.GAMMATOL ) THEN
         INFO = -12
      ELSE IF ( ABS( DELTA ).LT.DELTATOL ) THEN
         INFO = -13
      ELSE IF ( ABS( EPS ).LT.EPSTOL ) THEN
         INFO = -14
      ELSE IF ( ABS( XI ).LT.XITOL ) THEN
         INFO = -15
      ENDIF
*
*
      RLBL = -1
      IJOB = -1
*
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
*
      RETURN
*
*     End of QMRREVCOM
*
      END
