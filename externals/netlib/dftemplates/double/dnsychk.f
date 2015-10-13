*
      SUBROUTINE DNSYCHK( X, LDX, B, X0, WORK, LDW, PFORM, MATVEC,
     $                    MATVECTRANS, PSOLVE, PSOLVETRANS, PSOLVEQ,
     $                    PSOLVETRANSQ, BACKSOLVE, TOL, SCALEDTOL,
     $                    LTEST, NSYRES, NUMTESTS, NUMSUSP, CRITERR )
*
*     .. Scalar Arguments ..
      INTEGER            LDW, LDX, NUMTESTS, NUMSUSP, CRITERR
      DOUBLE PRECISION   TOL, SCALEDTOL
      LOGICAL            NSYRES
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( LDX,* ), B( * ), X0( * ), WORK( * )
      LOGICAL            LTEST( * )
      CHARACTER*5        PFORM( * )
*     ..
*     .. External Functions ..
      EXTERNAL         MATVEC, MATVECTRANS, PSOLVE, PSOLVETRANS,
     $                 PSOLVEQ, PSOLVETRANSQ, BACKSOLVE
*     ..
*     .. Common Blocks ..
*
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
*  Subroutine to test the performance of the nonstationary template
*  kernels on nonsymmetric matrices.
*
*  Generates, solves, and check accuracy of linear systems.
*
*  Algorithms tested:
*
*     4. BiCG
*     5. CGS
*     6. BiCGSTAB
*     7. GMRESm
*     8. QMR
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
*     .. Local Scalars ..
      INTEGER            I, J, MATFORM, NSYS, MAXIT, RESTRT, NPFORMS,
     $                   IPOINTR, IINDX, IVAL, NX, NY, NZ
      DOUBLE PRECISION   ANORM, MATNORM,
     $                   ONEFUN, ZEROFUN, NEGONEFUN, HENKFUN, HENKDFUN,
     $                   THOUSFUN, THOUSXFUN, TEN5X2FUN, NEGTHOUSXFUN
      LOGICAL            TSTBICG, TSTCGS, TSTBICGS, TSTGMRES, TSTQMR,
     $                   FLAG, LSAMEN
      CHARACTER          AFORM*4, RHS*4, INITIALGUESS*4
*     ..
*     .. Local Arrays ..
      INTEGER            INFO( 9 ), ITER( 9 )
      DOUBLE PRECISION   RESID( 9 )
*
*     PDE Coefficient functions.
*
      EXTERNAL           ONEFUN, ZEROFUN, NEGONEFUN, HENKFUN, HENKDFUN,
     $                   THOUSFUN, THOUSXFUN, TEN5X2FUN, NEGTHOUSXFUN
*
*     .. Executable Statements ..
*
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
      DO 5 I = 4, 8
         IF ( LTEST( I ) ) FLAG = .TRUE.
    5 CONTINUE
      IF ( .NOT.FLAG ) RETURN
*
      TSTBICG   = LTEST( 4 )
      TSTCGS    = LTEST( 5 )
      TSTBICGS  = LTEST( 6 )
      TSTGMRES  = LTEST( 7 )
      TSTQMR    = LTEST( 8 )
*
   10 CONTINUE
*
      DO 60 MATFORM = 1, NSYS
*
         READ(9,*) AFORM, NX, NY, NZ, RHS, INITIALGUESS
*
         N = NX*NY*NZ
         IPOINTR = 1
         IINDX   = IPOINTR + N+1
         IVAL    = IINDX   + 5*N
*
*        The following matrices are generated using a 5- or 7-point
*        stencil using centered differences on a 1d, 2d, or 3d grid,
*        with Dirichlet boundary conditions.
*
*        The last 7 arguments to this routine are the coefficient 
*        functions for the PDE:
*
*           delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
*           delx ( d u ) + dely (e u) + delz( f u ) + g u
*
         IF ( LSAMEN( 4, AFORM,'PDE1' ) ) THEN
*
*           u_xx + u_yy + au_x + (a_x/2)u for a = 20exp[3.5(x^2 + y^2 )]
*
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), ONEFUN, ONEFUN, ZEROFUN,
     $                    HENKFUN, ZEROFUN, ZEROFUN, HENKDFUN )
*
*        The following three PDE are from Yang, "PCG-like methods for
*        nonsymmetric linear systems".
*
         ELSE IF ( LSAMEN( 4, AFORM,'PDE2' ) ) THEN
*
*           u_xx + u_yy + u_zz + 1000u_x
*
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), ONEFUN, ONEFUN, ONEFUN,
     $                    THOUSFUN, ZEROFUN, ZEROFUN, ZEROFUN )
         ELSE IF ( LSAMEN( 4, AFORM,'PDE3' ) ) THEN
*
*           u_xx + u_yy + u_zz - 10^5x^2(u_x + u_y + u_z )
*
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), ONEFUN, ONEFUN, ONEFUN,
     $                    TEN5X2FUN, TEN5X2FUN, TEN5X2FUN, ZEROFUN )
         ELSE IF ( LSAMEN( 4, AFORM,'PDE4' ) ) THEN
*
*           u_xx + u_yy + u_zz + 1000exp(xyz)( u_x + u_y - u_z )
*
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), ONEFUN, ONEFUN, ONEFUN,
     $                    THOUSXFUN, THOUSXFUN, NEGTHOUSXFUN, ZEROFUN )
         ELSE
            WRITE(*,*) AFORM, 'IS AN UNKNOWM MATRIX TYPE'
            GO TO 60
         ENDIF
*
*        Convert to dense form.
*
         CALL COMP2DENSE( WORK(IVAL), WORK(IPOINTR), WORK(IINDX),
     $                    N, A, LDA,'ROW', INFO(1) )
         IF ( INFO(1).NE.0 ) THEN
            WRITE(*,81)  N, AFORM, INFO(1)
            WRITE(10,81) N, AFORM, INFO(1)
            GO TO 60
         ENDIF
*
   15    CONTINUE
*
         CALL VECGEN( RHS, N, A, LDA, B, INFO(2) )
         CALL VECGEN( INITIALGUESS, N, A, LDA, X0, INFO(3) )
*
         IF ( INFO(2).NE.0 ) THEN
            WRITE(*,92)  MATFORM
            WRITE(10,92) MATFORM
            GO TO 10
         ELSE IF ( INFO(3).NE.0 ) THEN
            WRITE(*,93)  MATFORM
            WRITE(10,93) MATFORM
            GO TO 10
         ENDIF
*
   20    CONTINUE
*
         MAXIT = 4*N
         ANORM = MATNORM( N, A, LDA )
*
*        Solve system using the various algorithms, using no
*        preconditioning, then diagonal preconditioning.
*
         DO 50 I = 1, NPFORMS
*
*           Initializations.
*
            DO 30 J = 3, 9
               INFO( J ) = 0
               ITER( J ) = MAXIT
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
               RESTRT = N / 2
               IF ( RESTRT.EQ.0 ) RESTRT = N
               CALL GMRES( N, B, X(1,7), RESTRT, WORK, LDW,
     $                     WORK((6+RESTRT)*LDW+1), LDW, ITER(7),
     $                     RESID(7), MATVEC, PSOLVE, INFO(7))
            ENDIF
*
            IF ( TSTQMR )
     $         CALL QMR( N, B, X(1,8), WORK, LDW, ITER(8), RESID(8),
     $                   MATVEC, MATVECTRANS, PSOLVEQ, PSOLVETRANSQ,
     $                   INFO(8) )
*
*           Check for convergence.
*
            DO 40 J = 4, 8
               IF ( INFO( J ).NE.0 ) NUMSUSP = NUMSUSP + 1
   40       CONTINUE
*
*           Write results to file.
*
            CALL RESULT( N, A, LDA, X, LDX, B, WORK(1), 'NSY', PFORM(I),
     $                   ITER, RESID, TOL, INFO, AFORM, ANORM,
     $                   LTEST, SCALEDTOL, NSYRES, CRITERR )
            NUMTESTS = NUMTESTS + 1
*
   50    CONTINUE
*
   60 CONTINUE
*
   81 FORMAT('WARNING: COULD NOT FORM ORDER ', I4, ' ', A4, 
     $       'MATRIX; INFO=',I2)
   92 FORMAT('ERROR: RHS', A4, ' NOT FORMED')
   93 FORMAT('ERROR: INITIAL GUESS', A4, ' NOT FORMED')
   94 FORMAT('ERROR: PRECONDITIONER', A5, ' NOT FORMED')
*
      RETURN
*
*     -- End of DNSYCHK
*
      END
