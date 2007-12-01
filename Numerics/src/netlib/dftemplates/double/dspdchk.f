*
      SUBROUTINE DSPDCHK( X, LDX, B, X0, WORK, LDW, PFORM, MATVEC, 
     $                    MATVECTRANS, PSOLVE, PSOLVETRANS, PSOLVEQ,
     $                    PSOLVETRANSQ, BACKSOLVE, TOL, SCALEDTOL,
     $                    LTEST, SPDRES, NUMTESTS, NUMSUSP, CRITERR )
*
*     .. Scalar Arguments ..
      INTEGER            LDW, LDX, NUMTESTS, NUMSUSP, CRITERR
      DOUBLE PRECISION   TOL, SCALEDTOL
      LOGICAL            SPDRES
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( LDX,* ), B( * ), X0( * ), WORK( * )
      LOGICAL            LTEST( * )
      CHARACTER*5        PFORM( * )
*     ..
*     .. External Functions ..
      EXTERNAL           MATVEC, MATVECTRANS, PSOLVE, PSOLVETRANS,
     $                   PSOLVEQ, PSOLVETRANSQ, BACKSOLVE
*     ..
*     .. Common Blocks ..
*     MAXDIM2 = MAXDIM*MAXDIM.
*
      INTEGER             MAXDIM, MAXDIM2
      PARAMETER         ( MAXDIM = 200, MAXDIM2 = 40000 )
*
      INTEGER             N, LDA
      DOUBLE PRECISION    A, M
      CHARACTER           CURPFORM*5
*
      COMMON            / SYSTEM / A( MAXDIM2 ), M( MAXDIM )
      COMMON            / MATDIM / N, LDA
      COMMON            / FORMS  / CURPFORM
*
*  Purpose
*  =======
*
*  Subroutine to test the performance of the template kernels
*  on symmetric positivie definite matrices.
*
*  Generates, solves, and checks accuracy of linear systems.
*
*  Algorithms tested:
*
*     1. CG
*     2. Cheby
*     3. SOR
*     4. BiCG
*     5. CGS
*     6. BiCGSTAB
*     7. GMRES
*     8. QMR
*     9. Jacobi ( for diagonally dominant Poisson matrices only )
*
*  Various systems are generated. Each method attempts to solve the 
*  system to the input TOL in MAXIT iterations. Each method iterates 
*  using various preconditioners. 
*
*  The result is compared with the normalized scaled residual 
*
*     || b - A*x || / ( ||A||||x||*N*TOL ). 
*
*  In order to do this, the solution vectors are stored in  matrix 
*  X( LDX,* ). Column j contains the solution vector for algorithm j,
*  j as defined above.
*
*  =================================================================
*
*     .. Parameters ..
*
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K, NX, NY, NZ, MATFORM, NSYS, NPFORMS,
     $                   MAXIT, RESTRT, IVAL, IINDX, IPOINTR
      DOUBLE PRECISION   ANORM, MATNORM,
     $                   ONEFUN, ZEROFUN, NEGONEFUN, HENKFUN, HENKDFUN,
     $                   THOUSFUN, THOUSXFUN, TEN5X2FUN, NEGTHOUSXFUN
      LOGICAL            TSTCG, TSTCHEBY, TSTSOR, TSTBICG, TSTCGS,
     $                   TSTBICGS, TSTGMRES, TSTQMR, TSTJACOBI,
     $                   FLAG, LSAMEN
      CHARACTER          AFORM*4, RHS*4, INITIALGUESS*4
*     ..
*     .. Local Arrays ..
      INTEGER            INFO( 9 ), ITER( 9 )
      DOUBLE PRECISION   RESID( 9 )
*     ..
*     .. External Functions ..
*     PDE Coefficient functions.
      EXTERNAL           ONEFUN, ZEROFUN, NEGONEFUN, HENKFUN, HENKDFUN,
     $                   THOUSFUN, THOUSXFUN, TEN5X2FUN, NEGTHOUSXFUN
*     ..
*     .. Executable Statements ..
*
      NPFORMS  = 2
      NUMTESTS = 0
      NUMSUSP  = 0
      CRITERR  = 0
*
      READ(9,*) NSYS
*
*     Check for quick return.
*
      IF ( NSYS.LT.0 ) THEN
         WRITE(*,*) 'ERROR IN NONSYMMETRIC TESTER: NUMBER OF SYSTEMS TO
     $BE GENERATED IS LESS THAN 0'
         RETURN
      ENDIF
*
      FLAG = .FALSE.
      DO 5 I = 1, 9
         IF ( LTEST( I ) ) FLAG = .TRUE.
    5 CONTINUE
      IF ( .NOT.FLAG ) RETURN
*
      TSTCG     = LTEST( 1 )
      TSTCHEBY  = LTEST( 2 )
      TSTSOR    = LTEST( 3 )
      TSTBICG   = LTEST( 4 )
      TSTCGS    = LTEST( 5 )
      TSTBICGS  = LTEST( 6 )
      TSTGMRES  = LTEST( 7 )
      TSTQMR    = LTEST( 8 )
      TSTJACOBI = LTEST( 9 )
*
   10 CONTINUE
*
      DO 60 MATFORM = 1, NSYS
*
         READ(9,*) AFORM, NX, NY, NZ, RHS, INITIALGUESS
*
*        The following two matrices are generated using a 5- or 7-point 
*        stencil using centered differences on a 1d, 2d, or 3d grid, 
*        with Dirichlet boundary conditions.
*
*        The last 7 arguments to this routine are the coefficient 
*        functions for the PDE:
*
*           delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
*           delx ( d u ) + dely (e u) + delz( f u ) + g u
*
         IF ( LSAMEN( 4, AFORM,'F2SH' ) ) THEN
*
*           -u_xx - u_yy
*
            N = NX*NY*NZ
            IPOINTR = 1
            IINDX   = IPOINTR + N+1
            IVAL    = IINDX   + 5*N
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), NEGONEFUN, NEGONEFUN,
     $                    ZEROFUN, ZEROFUN, ZEROFUN, ZEROFUN, ZEROFUN )
            CALL COMP2DENSE( WORK(IVAL), WORK(IPOINTR), 
     $                       WORK(IINDX), N, A, LDA,'ROW', INFO(1) )
            IF ( INFO(1).NE.0 ) THEN
               WRITE(*,81) N, INFO(1)
               WRITE(10,81) N, INFO(1)
               GO TO 60
            ENDIF
*
         ELSE IF ( LSAMEN( 4, AFORM,'F3SH' ) ) THEN
*
*           -u_xx - u_yy - u_zz
*
            N = NX*NY*NZ
            IPOINTR = 1
            IINDX   = IPOINTR + N+1
            IVAL    = IINDX   + 10*N
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), NEGONEFUN, NEGONEFUN,
     $                    NEGONEFUN, ZEROFUN, ZEROFUN, ZEROFUN,
     $                    ZEROFUN )
            CALL COMP2DENSE( WORK(IVAL), WORK(IPOINTR),
     $                       WORK(IINDX), N, A, LDA,'ROW', INFO(1) )
            IF ( INFO(1).NE.0 ) THEN
               WRITE(*,82) N, INFO(1)
               WRITE(10,82) N, INFO(1)
               GO TO 60
            ENDIF
         ELSE IF ( LSAMEN( 4, AFORM,'WATHEN' ) ) THEN
*
*           Form Wathen matrix. 
*           ( Matrix order will be 3*NX*NY + 2*NX + 2*NY + 1 )
*
            K = 0
            CALL WATHEN( NX, NY, K, N, A, LDA, WORK, LDW, INFO(1) )
            IF ( INFO(1).NE.0 ) THEN
               WRITE(*,83)  3*NX*NY + 2*NX + 2*NY + 1, INFO(1)
               WRITE(10,83) 3*NX*NY + 2*NX + 2*NY + 1, INFO(1)
               GO TO 60
            ENDIF
         ELSE
            WRITE(*,*) AFORM, 'IS AN UNKNOWM MATRIX TYPE'
            GO TO 60
         ENDIF
*
         MAXIT = 4*N
         ANORM = MATNORM( N, A, LDA )
*
*        Form RHS and initial guess vectors.
*
         CALL VECGEN( RHS, N, A, LDA, B, INFO(2) )
         CALL VECGEN( INITIALGUESS, N, A, LDA, X0, INFO(3) )
*
         IF ( INFO(2).NE.0 ) THEN
            WRITE(*,92)  RHS
            WRITE(10,92) RHS
            GO TO 60
         ELSE IF ( INFO(3).NE.0 ) THEN
            WRITE(*,93)  INITIALGUESS
            WRITE(10,93) INITIALGUESS
            GO TO 60
         ENDIF
*
   20    CONTINUE
*
*        Solve system using the various algorithms, using no 
*        preconditioning, then diagonal preconditioning.
*
         DO 50 I = 1, NPFORMS
*
            DO 30 J = 1, 9
               INFO( J ) = 0
               ITER( J )  = MAXIT
               RESID( J ) = TOL
               CALL DCOPY( N, X0, 1, X(1,J), 1 )
   30       CONTINUE
*
            CALL PRECON( N, A, LDA, PFORM( I ), M, INFO(1) )
            IF ( INFO(1).NE.0 ) THEN
               WRITE(*,94) PFORM(I)
               WRITE(10,94) PFORM(I)
               GO TO 50
            ENDIF
            CURPFORM = PFORM( I )
*
            IF ( TSTCG )
     $         CALL CG( N, B, X(1,1), WORK, LDW, ITER(1), RESID(1),
     $                  MATVEC, PSOLVE, INFO(1) )
*
            IF ( TSTCHEBY )
     $         CALL CHEBY( N, B, X(1,2), WORK, LDW, ITER(2), RESID(2),
     $                     MATVEC, PSOLVE, INFO(2) )
*
            IF  ( TSTSOR ) THEN
*
*              Set OMEGA
*
               IF ( MATFORM.EQ.1 ) THEN
*
*                 Guass-Seidel
*
                  WORK(1) = ONE
               ELSE
                  WORK(1) = 1.2D+0
               ENDIF
               CALL SOR( N, B, X(1,3), WORK, LDW, ITER(3),
     $                   RESID(3), MATVEC, BACKSOLVE, INFO(3) )
            ENDIF
*
            IF ( TSTBICG )
     $         CALL BICG( N, B, X(1,4), WORK, LDW, ITER(4), RESID(4),
     $                    MATVEC, MATVECTRANS, PSOLVE, PSOLVETRANS,
     $                    INFO(4) )
*
            IF ( TSTCGS )
     $         CALL CGS( N, B, X(1,5), WORK, LDW, ITER(5), RESID(5),
     $                   MATVEC, PSOLVE, INFO(5) )
*
            IF ( TSTBICGS )
     $         CALL BICGSTAB( N, B, X(1,6), WORK, LDW, ITER(6),
     $                        RESID(6), MATVEC, PSOLVE, INFO(6) )
*
            IF ( TSTGMRES ) THEN
*
*              For the symmetric case, restarts = N.
*
               RESTRT = N
               CALL GMRES( N, B, X(1,7), RESTRT, WORK, LDW,
     $                     WORK((6+RESTRT)*LDW+1), LDW, ITER(7),
     $                     RESID(7), MATVEC, PSOLVE, INFO(7) )
            ENDIF
*
            IF ( TSTQMR )
     $          CALL QMR( N, B, X(1,8), WORK, LDW, ITER(8), RESID(8),
     $                    MATVEC, MATVECTRANS, PSOLVEQ, PSOLVETRANSQ, 
     $                    INFO(8))
*
            IF ( TSTJACOBI ) THEN
               IF ( .NOT.LSAMEN( 3, AFORM,'WATH' ).AND.( I.EQ.1 ) ) THEN
*
*                 Since preconditioning does not apply to Jacobi, it is
*                 only called once per test matrix. (Wathen matrix is 
*                 not diagonally dominant.)
*
                  CALL JACOBI( N, B, X(1,9), WORK, LDW, ITER(9), 
     $                         RESID(9), MATVEC, INFO(9) )
                  IF ( INFO( 9 ).NE.0 ) NUMSUSP = NUMSUSP + 1
               ELSE
*
*                 Flag not to check accuracy.
*
                  INFO( 9 ) = 100
               ENDIF
            ENDIF
*
*           Check for convergence.
*
            DO 40 J = 1, 8
               IF ( INFO( J ).NE.0 ) NUMSUSP = NUMSUSP + 1 
   40       CONTINUE
*
*           Write results to file.
*
            CALL RESULT( N, A, LDA, X, LDX, B, WORK(1), 'SPD', PFORM(I),
     $                   ITER, RESID, TOL, INFO, AFORM, ANORM,
     $                   LTEST, SCALEDTOL, SPDRES, CRITERR )
            NUMTESTS = NUMTESTS + 1
*
   50    CONTINUE
*
   60 CONTINUE
*
   81 FORMAT('WARNING: COULD NOT FORM 2-D POISSON (N=', I4,'),INFO=',I2)
   82 FORMAT('WARNING: COULD NOT FORM 3-D POISSON (N=', I4,'),INFO=',I2)
   83 FORMAT('WARNING: COULD NOT FORM ORDER', I4,' WATHEN MATRIX')
*
   92 FORMAT('ERROR: RHS FOR TEST MATRIX', A4, ' NOT FORMED')
   93 FORMAT('ERROR: INITIAL GUESS FOR TEST MATRIX', A4, ' NOT FORMED')
   94 FORMAT('ERROR: PRECONDITIONER', A5, ' NOT FORMED')
*
  200 CONTINUE
*
      RETURN
*
*     -- End of DSPDCHK
*
      END
