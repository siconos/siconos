*
      PROGRAM TEMPLATESTESTER
*
*  Test program for the DOUBLE PRECISION iterative templates.
*
*  The program must be driven by a short data file. The first 18 records
*  of the file are read using list-directed input, the last 16 records
*  are read using the format ( A6, L2 ). An annotated example of a data
*  file s a follows:
*
*  1.0D-15                        CONVERGENCE TOLERANCE
*  10                             SCALED RESIDUAL TOLERANCE
*  CG     T PUT F FOR NO TEST.    ALGORITHMS TO BE TESTED
*  CHEBY  T PUT F FOR NO TEST.
*  SOR    T PUT F FOR NO TEST.
*  BICG   T PUT F FOR NO TEST.
*  CGS    T PUT F FOR NO TEST.
*  BICGS  T PUT F FOR NO TEST.
*  GMRES  F PUT F FOR NO TEST.
*  QMR    T PUT F FOR NO TEST.
*  JACOB  T PUT F FOR NO TEST.
*  3                              NUMBER OF SPD MATRICES TO BE GENERATED
*  WATH  2, 2, 1, ONES, ZERO      MATRIX, NX, NY, NZ, RHS, INITIAL GUESS
*  F2SH  6, 6, 1, SUMR, ZERO
*  F3SH  3, 3, 3, ONES, ZERO
*  BICG   T  PUT F FOR NO TEST.   ALGORITHMS TO BE TESTED
*  CGS    T  PUT F FOR NO TEST.
*  BICGS  T  PUT F FOR NO TEST.
*  GMRES  F  PUT F FOR NO TEST.
*  QMR    T  PUT F FOR NO TEST.
*  4                              NUMBER OF MATTRICES TO BE GENERATED
*  PDE1, 5, 5, 5, SUMR , ZERO     MATRIX, NX, NY, NZ, RHS, INITIAL GUESS
*  PDE2, 5, 5, 5, SUMR , ZERO
*  PDE3, 5, 5, 5, ONES , ZERO
*  PDE4, 6, 6, 1, ONES , ZERO
*
*  See:
*
*     Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     Eijkhout, Pozo, Romine, and van der Vorst.
*     Templates for the Solution of Linear Systems: Building Blocks
*     for Iterative Methods, SIAM Publications, 1993. 
*     (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*  -- Written on 1-November-1993.
*     Richard Barrett, University of Tennessee
*     Jack Dongarra, Univ. of Tennessee and Oak Ridge National 
*     Laboratory
*
*     .. Parameters ..
*     MAXLEN must be greater than or equal to (2N**2)+3N, i.e. WORK
*     must have dimension N x (2N+3). This is for SOR (see StatUtils
*     for details). Chebyshev requires N*2. For workspace requirements
*     of the algorithms, see the individial template.
*
      INTEGER            MAXDIM, MAXLEN, NSUBS
      PARAMETER        ( MAXDIM = 200, MAXLEN = 80000, NSUBS = 9)
*     ..
*     .. Scalar Declarations ..
      INTEGER            I, LDX, LDW, SPDTESTS, NSYTESTS, SUSPSPD,
     $                   SUSPNSY, CRITSPD, CRITNSY
      DOUBLE PRECISION   TOL, SCALEDTOL
      LOGICAL            LTESTT, LSAMEN, SPDRES, NSYRES
      CHARACTER*6        SNAMET
*     ..
*     .. Array Declarations ..
      DOUBLE PRECISION   X( MAXDIM,NSUBS ), B( MAXDIM ), X0( MAXDIM ),
     $                   WORK( MAXLEN )
      LOGICAL            LTEST( NSUBS )
      CHARACTER*6        SNAMES( NSUBS )
      CHARACTER*5        PFORM( 2 )
*     ..
*     .. Common Blocks ..
      INTEGER            N, LDA
      COMMON           / MATDIM / N, LDA
*     ..
*     .. External Routines ..
      EXTERNAL           MATVEC, MATVECTRANS, PSOLVE, PSOLVETRANS,
     $                   PSOLVEQ, PSOLVETRANSQ, BACKSOLVE
*
      DATA               SNAMES/'CG    ', 'CHEBY ', 'SOR   ', 'BICG  ',
     $                   'CGS   ', 'BICGS ', 'GMRES ', 'QMR   ',
     $                   'JACOB '/
*     ..
*     .. Executable Statements ..
*
*     Initializations.
*
      LDA = MAXDIM
      LDX = MAXDIM
      LDW = MAXDIM
*
      SPDRES = .TRUE.
      NSYRES = .TRUE.
*
      PFORM( 1 ) = 'IDENT'
      PFORM( 2 ) = 'JACBI'
*
      OPEN( UNIT = 9, FILE = 'test.data' )
      OPEN( UNIT = 10, FILE = 'test.results' )
*
*     Get the convergence tolerance, the tolerance for the normalized
*     scaled residual, and the number of systems to be generated.
*     and the algorithms to be tested.
*
      READ(9,*) TOL
      READ(9,*) SCALEDTOL
*
*     Get input data for SPD testing: 
*     Read names of subroutines and flags which indicate whether 
*     they are to be tested.
*
      DO 10 I = 1, NSUBS
         LTEST( I ) = .FALSE.
   10 CONTINUE
   20 READ( 9, FMT = 998 )SNAMET, LTESTT
      DO 30 I = 1, NSUBS
         IF( LSAMEN( 6, SNAMET, SNAMES( I ) ) ) GO TO 40
   30 CONTINUE
      WRITE( *, FMT = 999 )SNAMET
      STOP
   40 LTEST( I ) = LTESTT
      IF ( I.LT.NSUBS ) GO TO 20
*
   50 CONTINUE
*
*     Begin testing.
*
      CALL HEADER( TOL )
*
*     Symmetric Positive Definite Routine Tester.
*
      CALL DSPDCHK( X, LDX, B, X0, WORK, LDW, PFORM, MATVEC,
     $              MATVECTRANS, PSOLVE, PSOLVETRANS, PSOLVEQ,
     $              PSOLVETRANSQ, BACKSOLVE, TOL, SCALEDTOL, LTEST,
     $              SPDRES, SPDTESTS, SUSPSPD, CRITSPD )
*
*     Get input data for Nonsymmetric testing:
*     Read names of subroutines and flags which indicate whether 
*     they are to be tested.
*     
      DO 60 I = 1, NSUBS
         LTEST( I ) = .FALSE.
   60 CONTINUE
   70 READ( 9, FMT = 998, END = 100 )SNAMET, LTESTT
      DO 80 I = 4, 8
         IF( LSAMEN( 6, SNAMET, SNAMES( I ) ) ) GO TO 90
   80 CONTINUE
      WRITE( *, FMT = 999 )SNAMET
      STOP
   90 LTEST( I ) = LTESTT
      IF ( I.LT.8 ) GO TO 70
*
  100 CONTINUE
*
*     Nonsymmetric Routine Tester.
*
      CALL DNSYCHK( X, LDX, B, X0, WORK, LDW, PFORM, MATVEC,
     $              MATVECTRANS, PSOLVE, PSOLVETRANS, PSOLVEQ,
     $              PSOLVETRANSQ, BACKSOLVE, TOL, SCALEDTOL, LTEST, 
     $              NSYRES, NSYTESTS, SUSPNSY, CRITNSY )
*
*     End of testing.
*
      CALL FOOTER()
*
      CLOSE( UNIT =  9 )
      CLOSE( UNIT = 10 )
*
*     Print overall results to screen.
*
      WRITE(*,*)
      IF ( ( SPDRES ).AND.( NSYRES ) ) THEN
*
*        All tests passed.
*
         WRITE(*,*) 'TESTS COMPLETE:'
         WRITE(*,*)
         IF ( SPDTESTS.GT.0 ) WRITE(*,900) SPDTESTS
         IF ( NSYTESTS.GT.0 ) WRITE(*,901) NSYTESTS
      ELSE 
         IF ( SPDRES ) THEN
*
            IF ( SPDTESTS.GT.0 ) THEN
*
*              SPD tests passed.
*
               WRITE(*,*) 'TESTS COMPLETE:'
               WRITE(*,*)
               WRITE(*,910) SPDTESTS
            ELSE
               WRITE(*,*)
               WRITE(*,*) 'SPD TESTING NOT PERFORMED.'
            ENDIF
         ELSE
*
*           SPD testing failed.
*
            WRITE(*,911) SPDTESTS
            WRITE(*,990) CRITSPD, SUSPSPD
            WRITE(*,991)
         ENDIF
         WRITE(*,*)
         IF ( NSYRES ) THEN
*
            IF ( SPDTESTS.GT.0 ) THEN
*
*              Nonsymmetric tests passed.
*
               WRITE(*,*) 'TESTS COMPLETE:'
               WRITE(*,*)
               WRITE(*,920) NSYTESTS
            ELSE
               WRITE(*,*)
               WRITE(*,*) 'NONSYMMETRIC TESTING NOT PERFORMED.'
            ENDIF
         ELSE
*
*           Nonsymmetric testing failed.
*
            WRITE(*,*)
            WRITE(*,921) NSYTESTS
            WRITE(*,990) CRITNSY, SUSPNSY
            WRITE(*,991)
         ENDIF
      ENDIF
      WRITE(*,*)
*
*     Format statements for screen output of general test results.
*
  900 FORMAT('   SYMMETRIC POSITIVE DEFINITE ROUTINES PASSED. (', I3,
     $       ' TESTS )')
  901 FORMAT('   NONSYMMETRIC ROUTINES PASSED. (', I3,' TESTS )')
*
  910 FORMAT(' PASSED FOR SYMMETRIC POSITIVE DEFINITE MATRICES (' ,I3,' 
     $TESTS )')
  911 FORMAT(' SYMMETRIC POSITIVE DEFINITE MATRICES: (' , I3,' TESTS )')
*
  920 FORMAT(' PASSED FOR NONSYMMETRIC MATRICES (' ,I3,' TESTS )')
  921 FORMAT(' NONSYMMETRIC MATRICES: (' , I3,' TESTS )')
*
  990 FORMAT(' THERE ARE',I3,' CRITICAL ERRORS AND ' ,I3,' SUSPICIOUS RE
     $SULTS.')
  991 FORMAT(' SEE FILE test.results FOR ADDITIONAL INFORMATION')
*
  998 FORMAT( A6, L2 )
  999 FORMAT( ' SUBPROGRAM NAME ', A6, ' NOT RECOGNIZED', /' ******* T',
     $      'ESTS ABANDONED *******' )
*
      STOP
*
*     End of Driver for Testing the Iterative Templates
*
      END
*
*     ===============================================================
      SUBROUTINE VECGEN( FORM, N, A, LDA, B, INFO )
*
      DOUBLE PRECISION           ZERO, ONE
      PARAMETER    ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*
      INTEGER        I, J, N, LDA, INFO
      DOUBLE PRECISION           A( LDA,* ), B( * ), TMP
      CHARACTER      FORM*4
      LOGICAL        LSAMEN
*
      INFO = 0
*
      IF ( LSAMEN( 3, FORM,'ONES' ) ) THEN
         DO 10 I = 1, N
            B( I ) = ONE
   10    CONTINUE
      ELSE IF ( LSAMEN( 3, FORM,'ZEROS' ) ) THEN
         DO 20 I = 1, N
            B( I ) = ZERO
   20    CONTINUE
      ELSE IF ( LSAMEN( 3, FORM,'SUMROW' ) ) THEN
         DO 40 I = 1, N
            TMP = ZERO
            DO 30 J = 1, N
               TMP = TMP + A( I,J )
   30       CONTINUE
            B( I ) = TMP
   40    CONTINUE
      ELSE
         INFO = -1
      ENDIF
*
      RETURN
*
      END
*
*     ===============================================================
      SUBROUTINE PRECON( N, A, LDA, PFORM, M, INFO )
*
*     .. Scalar and Array Declarations ..
*
      INTEGER       N, LDA, INFO
      DOUBLE PRECISION          A( LDA,* ), M( * )
      CHARACTER *5  PFORM
*
*  Purpose:
*  =======
*
*  PRECON forms a preconditioner matrix of type PROFRM for 
*  iterative solvers of the linear system Ax = b.
*
*  PFORM:
*
*      IDENT        identity matrix (for testing)
*
*      JACBI        diagonal scaling
*
*     ==============================================
*
*     .. Local Scalars ..
*
      INTEGER             I
      LOGICAL             LSAMEN
*
*     .. Executable Statements ..
*
      IF ( LSAMEN( 5, PFORM,'IDENT' ) ) THEN
*
*         Identity matrix need not be formed, since the solve involving
*         the preconditioner (PSolve) merely copies the right hand side
*         to the solution vector.
*
          RETURN
*
      ELSE IF ( LSAMEN( 5, PFORM,'JACBI' ) ) THEN
*
*         Diagonal Scaling: diag(A). Note that we actually form inv(M) so that
*         solver can use multiplication.
*
          DO 10 I = 1, N
             M( I ) = A( I,I )
   10     CONTINUE
*
      ELSE
*
*        Selected preconditioner not implemented
*
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*) 'PRECONDITIONER ',PFORM,' NOT YET IMPLEMENTED'
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)
         INFO = -1
*
      ENDIF
*
      RETURN
      END
* 
*     ================================================================
      DOUBLE PRECISION FUNCTION GETBREAK()
*
*     Get breakdown parameter tolerance; for the test routine,
*     set to machine precision.
*
      DOUBLE PRECISION EPS, DLAMCH
*
      EPS = DLAMCH('EPS')
      GETBREAK = EPS**2
*
      RETURN
*
      END
*     ===============================================================
      DOUBLE PRECISION FUNCTION SCALEDRESID( ANORM, N, X, RK, TOL )
*
*     Returns |B-A*X| / ( |A||X|*N*TOL ), using the infinity norm.
*
      INTEGER            N, IDAMAX
      DOUBLE PRECISION   ANORM, TOL, XNORM, RESNORM,
     $                   X( * ), RK( * )
*
      XNORM   = ABS( X( IDAMAX( N, X, 1 ) ) )
      RESNORM = ABS( RK( IDAMAX( N, RK, 1 ) ) )
*
      SCALEDRESID = RESNORM / ( TOL * N * ANORM * XNORM )
*
      RETURN
*
      END
*     ===========================================================
      DOUBLE PRECISION FUNCTION MATNORM( N, A, LDA )
*
*     Compute infinity norm of matrix A.
*
      INTEGER            N, LDA, I, J
      DOUBLE PRECISION   ROWSUM, ZERO, TEMP, A( LDA,* )
      PARAMETER        ( ZERO = 0.0D+0 )
*
      TEMP = ZERO
      DO 20 I = 1, N
         ROWSUM = ZERO
         DO 10 J = 1, N
            ROWSUM = ROWSUM + ABS( A( I,J ) )
   10    CONTINUE
         TEMP = MAX( ROWSUM, TEMP )
   20 CONTINUE
*
      MATNORM = TEMP
*
      RETURN
*
      END
*
*     ===============================================================
      SUBROUTINE RESULT( N, A, LDA, X, LDX, B, RK, MATTYPE, PFORM,
     $                   ITER, RESID, TOL, INFO, AFORM, ANORM, 
     $                   LTEST, SCALEDTOL, TESTPASSED, CRITERR )
*
*     .. Argument Declaractions ..
*
      INTEGER             N, LDA, LDX, CRITERR,
     $                    ITER( * ), INFO( * )
      DOUBLE PRECISION    TOL, ANORM, SCALEDTOL,
     $                    A( LDA,* ), X( LDX,* ), B( * ), RK( * ),
     $                    RESID( * )
      CHARACTER           MATTYPE*3, AFORM*4, PFORM*5
      LOGICAL             TESTPASSED, LTEST( * )
*
*  Purpose
*  =======
*
*  Report results of METHOD on matrix type MATTYPE. If the residual 
*  is not directly computed by the algorithm, then the residual RESID 
*  as returned by the algorithm is compared with the residual as 
*  computed using the solution returned, i.e. || B-AX ||. 
*  =======================================================
*
*     .. Local Declarations ..
*
      DOUBLE PRECISION    ZERO, ONE
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*
      INTEGER             I, FIRSTALG, NUMALG
      DOUBLE PRECISION    DNRM2, SCALEDRESID, SRESID( 9 )
      CHARACTER*9         METHOD( 9 )
      LOGICAL             LSAME, LSAMEN
*
      EXTERNAL            DCOPY, DGEMV, DNRM2
*
*     .. Executable Statements ..
*
      METHOD( 1 ) = 'CG       '
      METHOD( 2 ) = 'Chebyshev'
      METHOD( 3 ) = 'SOR      '
      METHOD( 4 ) = 'BiCG     '
      METHOD( 5 ) = 'CGS      '
      METHOD( 6 ) = 'BiCGSTAB '
      METHOD( 7 ) = 'GMRESm   '
      METHOD( 8 ) = 'QMR      '
      METHOD( 9 ) = 'Jacobi   '
*
*     Compare algorithm reported residual with |b-AX|/(|A||x|n*TOL)
*
      IF ( LSAME( MATTYPE,'SPD' ) ) THEN
         FIRSTALG = 1
         NUMALG   = 9
      ELSE
         FIRSTALG = 4
         NUMALG   = 8
      ENDIF 
      DO 10 I = FIRSTALG, NUMALG
         IF ( RESID( I ).NE.ZERO ) THEN
            CALL DCOPY( N, B, 1, RK, 1 )
            CALL DGEMV('N', N, N, -ONE, A, LDA, X( 1,I ), 1, ONE,
     $                  RK, 1 )
            SRESID( I ) = SCALEDRESID( ANORM, N, X( 1,I ), RK, TOL )
         ENDIF
   10 CONTINUE
*
      IF ( LSAMEN( 4, AFORM,'F2SH') ) THEN
         IF ( LSAMEN( 2, PFORM,'IDENT' ) ) THEN
            WRITE(10,900) N
            WRITE(*,900) N
         ELSE IF ( LSAMEN( 2, PFORM,'JACBI' ) ) THEN
            WRITE(10,901) N
            WRITE(*,901) N
         ENDIF
      ELSE IF ( LSAMEN( 4, AFORM,'F3SH') ) THEN
         IF ( LSAMEN( 2, PFORM,'IDENT' ) ) THEN
            WRITE(10,902) N
            WRITE(*,902) N
         ELSE IF ( LSAMEN( 2, PFORM,'JACBI' ) ) THEN
            WRITE(10,903) N
            WRITE(*,903) N
         ENDIF
      ELSE IF ( LSAMEN( 4, AFORM,'WATH' ) ) THEN
         IF ( LSAMEN( 2, PFORM,'IDENT' ) ) THEN
            WRITE(10,910) N
            WRITE(*,910) N
         ELSE IF ( LSAMEN( 2, PFORM,'JACBI' ) ) THEN
            WRITE(10,911) N
            WRITE(*,911) N
         ENDIF
      ELSE IF ( LSAMEN( 3, AFORM,'PDE' ) ) THEN
         IF ( LSAMEN( 2, PFORM,'IDENT' ) ) THEN
            WRITE(10,920) N, AFORM
            WRITE(*,920) N, AFORM
         ELSE IF ( LSAMEN( 2, PFORM,'JACBI' ) ) THEN
            WRITE(10,921) N, AFORM
            WRITE(*,921) N, AFORM
         ENDIF
      ENDIF
      WRITE(10,*)
*
*     Loop over the algorithms, with a final error check.
*
      DO 30 I = FIRSTALG, NUMALG
*
*        Check updated residual vs. scaled residual.
*
         IF ( LTEST( I ) ) THEN
            IF ( INFO( I ).EQ.0 ) THEN
*
*              Method claims to have found solution. 
*              Check scaled residual.
*
               IF ( SRESID( I ).LE.SCALEDTOL ) THEN
*
*                 Scaled residual check passed.
*
                  WRITE(10,991) METHOD( I ), RESID( I ), SRESID( I ),
     $                          ITER( I )
                  WRITE(*,991) METHOD( I ), RESID( I ), SRESID( I ),
     $                         ITER( I )
               ELSE
                  CRITERR = CRITERR + 1
                  TESTPASSED = .FALSE.
                  WRITE(10,992) METHOD( I ), RESID( I ), SRESID( I ),
     $                          ITER( I )
                  WRITE(*,992) METHOD( I ), RESID( I ), SRESID( I ),
     $                         ITER( I )
               ENDIF
            ELSE IF ( INFO( I ).EQ.100 ) THEN
               GO TO 30
            ELSE
               TESTPASSED = .FALSE.
*
*              Method claims to have not found solution to tolerance,
*              either because the maximum number of iterations were
*              performed, or breakdown occured.
*
               WRITE(10,993) METHOD( I ), RESID( I ), SRESID( I ),
     $                       ITER( I ), INFO( I )
               WRITE(*,993) METHOD( I ), RESID( I ), SRESID( I ),
     $                       ITER( I ), INFO( I )
            ENDIF
*
         ELSE
*
*           Method was not involved in test
*
            GO TO 30 
*
         ENDIF
*
   30 CONTINUE
*
      WRITE(10,*) '-----------------------------------------------------
     $--'
*
*     Header for each system.
*
  900 FORMAT('Order', I4,' SPD 2-d Poisson matrix (no preconditioning)')
  901 FORMAT('Order', I4,' SPD 2-d Poisson matrix (Jacobi preconditionin
     $g)')
  902 FORMAT('Order', I4,' SPD 3-d Poisson matrix (no preconditioning)')
  903 FORMAT('Order', I4,' SPD 3-d Poisson matrix (Jacobi preconditionin
     $g)')
  910 FORMAT('Order ', I4,' SPD Wathen matrix (no preconditioning)')
  911 FORMAT('Order ', I4,' SPD Wathen matrix (Jacobi preconditioning)')
  920 FORMAT('Order ', I4,' ', A4, ' nonsymmetric matrix (no preconditio
     $ning)')
  921 FORMAT('Order ', I4,' ', A4, ' nonsymmetric matrix (Jacobi precond
     $itioning)')
*
*     Reporting of results.
*
  991 FORMAT('  ', A9,' ',1PE8.2,'    ', 1PE8.2,'   ', I5 )
  992 FORMAT('  ', A9,' ',1PE8.2,'    ', 1PE8.2,'   ', I5, 
     $        '            X' )
  993 FORMAT('  ', A9,' ',1PE8.2,'    ', 1PE8.2,'   ', I5,'    ', I3 )
*
      RETURN
*
*     End of Result.f
*
      END
*
*     ==================================================================
      SUBROUTINE HEADER( TOL )
*
      DOUBLE PRECISION   TOL, EPS, DLAMCH
*
      EPS = DLAMCH('E')
*
*     Print header to file.
*
      WRITE(10,*)
      WRITE(10,*) 'DETAILS OF ITERATIVE TEMPLATES TEST:'
      WRITE(10,*)
      WRITE(10,*) '   Univ. of Tennessee and Oak Ridge National Laborato
     $ry'
      WRITE(10,*) '   October 1, 1993'
      WRITE(10,*) '   Details of these algorithms are described in "Temp
     $lates'
      WRITE(10,*) '   for the Solution of Linear Systems: Building Block
     $s for'
      WRITE(10,*) '   Iterative Methods", Barrett, Berry, Chan, Demmel, 
     $Donato,'
      WRITE(10,*) '   Dongarra, Eijkhout, Pozo, Romine, and van der Vors
     $t,'
      WRITE(10,*) '   SIAM Publications, 1993.'
      WRITE(10,*) '   (ftp netlib2.cs.utk.edu; cd linalg; get templates.
     $ps).'
      WRITE(10,*)
      WRITE(10,*)
      WRITE(10,21) EPS
      WRITE(10,22) TOL
      WRITE(10,*)
      WRITE(10,*)
      WRITE(10,*) ' For a detailed description of the following informat
     $ion,'
      WRITE(10,*) ' see the end of this file.'
      WRITE(10,*)
      WRITE(10,*) '=====================================================
     $='
      WRITE(10,*) '           CONVERGENCE  NORMALIZED  NUM'
      WRITE(10,*) '  METHOD    CRITERION    RESIDUAL   ITER  INFO  FLAG'
      WRITE(10,*) '=====================================================
     $='
      WRITE(10,*)
*
*     Print header to screen.
*
      WRITE(*,*)
      WRITE(*,*) 'DETAILS OF ITERATIVE TEMPLATES TEST:'
      WRITE(*,*)
      WRITE(*,*) '   Univ. of Tennessee and Oak Ridge National Laborator
     $y'
      WRITE(*,*) '   October 1, 1993'
      WRITE(*,*) '   Details of these algorithms are described in "Templ
     $ates'
      WRITE(*,*) '   for the Solution of Linear Systems: Building Blocks
     $ for'
      WRITE(*,*) '   Iterative Methods", Barrett, Berry, Chan, Demmel, D
     $onato,'
      WRITE(*,*) '   Dongarra, Eijkhout, Pozo, Romine, and van der Vorst
     $,'
      WRITE(*,*) '   SIAM Publications, 1993.'
      WRITE(*,*) '   (ftp netlib2.cs.utk.edu; cd linalg; get templates.p
     $s).'
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,21) EPS
      WRITE(*,22) TOL
      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) ' For a detailed description of the following informati
     $on,'
      WRITE(10,*) ' see the end of this file.'
      WRITE(*,*)
      WRITE(*,*) '====================================================='
      WRITE(*,*) '           CONVERGENCE  NORMALIZED  NUM'
      WRITE(*,*) '  METHOD    CRITERION    RESIDUAL   ITER  INFO FLAG'
      WRITE(*,*) '====================================================='
      WRITE(*,*)
*
  21  FORMAT( 'MACHINE PRECISION = ', 1PE8.2 )
  22  FORMAT( 'CONVERGENCE TEST TOLERANCE = ', 1PE8.2 )

*
      RETURN
*
      END
*
*     ==================================================================
      SUBROUTINE FOOTER()
*
*     Puts descriptive information at bottom of results file
*
      WRITE(10,*)
      WRITE(10,*) '======'
      WRITE(10,*) 'LEGEND:'
      WRITE(10,*) '======'
      WRITE(10,*)
      WRITE(10,*) '   =================='
      WRITE(10,*) '   SYSTEM DESCRIPTION'
      WRITE(10,*) '   =================='
      WRITE(10,*)
      WRITE(10,*) '   SPD matrices:'
      WRITE(10,*)
      WRITE(10,*) '      WATH: "Wathen Matrix": consistent mass matrix'
      WRITE(10,*) '      F2SH: 2-d Poisson problem'
      WRITE(10,*) '      F3SH: 3-d Poisson problem'
      WRITE(10,*)
      WRITE(10,*) '      PDE1: u_xx+u_yy+au_x+(a_x/2)u'
      WRITE(10,*) '            for a = 20exp[3.5(x**2+y**2 )]'
      WRITE(10,*)
      WRITE(10,*) '   Nonsymmetric matrices:'
      WRITE(10,*)
      WRITE(10,*) '      PDE2: u_xx+u_yy+u_zz+1000u_x'
      WRITE(10,*) '      PDE3  u_xx+u_yy+u_zz-10**5x**2(u_x+u_y+u_z )'
      WRITE(10,*) '      PDE4: u_xx+u_yy+u_zz+1000exp(xyz)(u_x+u_y-u_z)'
      WRITE(10,*)
      WRITE(10,*) '   ====================='
      WRITE(10,*) '   CONVERGENCE CRITERION'
      WRITE(10,*) '   ====================='
      WRITE(10,*)
      WRITE(10,*) '   Convergence criteria: residual as reported by the'
      WRITE(10,*) '   algorithm: ||AX - B|| / ||B||. Note that NaN may s
     $ignify'
      WRITE(10,*) '   divergence of the residual to the point of numeric
     $al overflow.'
      WRITE(10,*) 
      WRITE(10,*) '   ==================='
      WRITE(10,*) '   NORMALIZED RESIDUAL'
      WRITE(10,*) '   ==================='
      WRITE(10,*)
      WRITE(10,*) '   Normalized Residual: ||AX - B|| / (||A||||X||*N*TO
     $L).'
      WRITE(10,*) '   This is an apostiori check of the iterated solutio
     $n.'
      WRITE(10,*)
      WRITE(10,*) '   ===='
      WRITE(10,*) '   INFO'
      WRITE(10,*) '   ===='
      WRITE(10,*) 
      WRITE(10,*) '   If this column is blank, then the algorithm claims
     $ to have'
      WRITE(10,*) '   found the solution to tolerance (i.e. INFO = 0).'
      WRITE(10,*) '   This should be verified by checking the normalized
     $residual.'
      WRITE(10,*)
      WRITE(10,*) '   Otherwise:'
      WRITE(10,*)
      WRITE(10,*) '      = 1: Convergence not achieved given the maximum
     $ number of iterations.'
      WRITE(10,*)
      WRITE(10,*) '      Input parameter errors:'
      WRITE(10,*) 
      WRITE(10,*) '      = -1: matrix dimension N < 0'
      WRITE(10,*) '      = -2: LDW < N'
      WRITE(10,*) '      = -3: Maximum number of iterations <= 0.'
      WRITE(10,*) '      = -4: For SOR: OMEGA not in interval (0,2)'
      WRITE(10,*) '            For GMRES: LDW2 < 2*RESTRT'
      WRITE(10,*) '      = -5: incorrect index request by uper level.'
      WRITE(10,*) '      = -6: incorrect job code from upper level.'
      WRITE(10,*)
      WRITE(10,*) '      <= -10: Algorithm was terminated due to breakdo
     $wn.'
      WRITE(10,*) '              See algorithm documentation for details
     $.'
      WRITE(10,*)
      WRITE(10,*) '   ===='
      WRITE(10,*) '   FLAG'
      WRITE(10,*) '   ===='
      WRITE(10,*)
      WRITE(10,*) '      X: Algorithm has reported convergence, but'
      WRITE(10,*) '         approximate solution fails scaled'
      WRITE(10,*) '         residual check.'
      WRITE(10,*)
      WRITE(10,*) '   ====='
      WRITE(10,*) '   NOTES'
      WRITE(10,*) '   ====='
      WRITE(10,*)
      WRITE(10,*) '   GMRES: For the symmetric test matrices, the restar
     $t parameter is'
      WRITE(10,*) '   set to N. This should, theoretically, result in no
     $ restarting. For'
      WRITE(10,*) '   nonsymmetric testing the restart parameter is set 
     $to N / 2.'
      WRITE(10,*)
      WRITE(10,*) '   Stationary methods:'
      WRITE(10,*)
      WRITE(10,*) '   - Since the residual norm ||b-Ax|| is not availabl
     $e as part of'
      WRITE(10,*) '     the algorithm, the convergence criteria is diffe
     $rent from the'
      WRITE(10,*) '     nonstationary methods. Here we use'
      WRITE(10,*)
      WRITE(10,*) '        || X - X1 || / || X ||.'
      WRITE(10,*) 
      WRITE(10,*) '     That is, we compare the current approximated sol
     $ution with the'
      WRITE(10,*) '     approximation from the previous step.'
      WRITE(10,*)
      WRITE(10,*) '   - Since Jacobi and SOR do not use preconditioning,
     $'
      WRITE(10,*) '     Jacobi is only iterated once per system, and SOR
     $ loops over'
      WRITE(10,*) '     different values for OMEGA (the first time throu
     $gh OMEGA = 1,'
      WRITE(10,*) '     i.e. the algorithm defaults to Gauss-Siedel). Th
     $is explains the '
      WRITE(10,*) '     different residual norms for SOR with the same m
     $atrix.'
      WRITE(10,*)
*
      RETURN
*
      END
