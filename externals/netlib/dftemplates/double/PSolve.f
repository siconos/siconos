*
*  This file contains the preconditioner solve routines:
*
*  PSOLVE and PSOLVETRANS call the appropriate solver:
*
*     PSOLVENONE and PSOLVENONETRANS for using no preconditioning.
*
*     PSOLVEJAC and PSOLVEJACTRANS for Jacobi preconditioning.
*
*     Also included are the solvers for QMR which require left and right
*     preconditioning: PSOLVEQ and PSOLVETRANSQ
*
*     ================================================================
      SUBROUTINE PSOLVE( X, B )
*
*     .. Array Arguments ..
      DOUBLE PRECISION    X( * ), B( * )
*     ..
*     .. Local Scalars ..
      LOGICAL             LSAME
*     ..
*     .. Common Blocks ..
      CHARACTER           CURPFORM*5
      COMMON            / FORMS  / CURPFORM
*
*     .. Executable Statements ..
*
      IF ( LSAME( CURPFORM, 'IDENT' ) ) THEN
         CALL PSOLVENONE( X, B )
      ELSE IF ( LSAME( CURPFORM, 'JACBI' ) ) THEN
         CALL PSOLVEJAC( X, B )
      ELSE
         WRITE(*,*) 'IN PSOLVE: UNKNOWN PRECONDITIONER', CURPFORM, 
     $              ' QUITTING' 
         STOP
      ENDIF
*
      RETURN
*
      END
*     ================================================================
*
      SUBROUTINE PSOLVETRANS( X, B )
*
*     .. Array Arguments ..
      DOUBLE PRECISION    X( * ), B( * )
*     ..
*     .. Common Blocks ..
      CHARACTER           CURPFORM*5
      COMMON            / FORMS  / CURPFORM
*     ..
*     .. Local Scalars ..
      LOGICAL          LSAME
*
*     .. Executable Statements ..
*
      IF ( LSAME( CURPFORM, 'IDENT' ) ) THEN
         CALL PSOLVENONETRANS( X, B )
      ELSE IF ( LSAME( CURPFORM, 'JACBI' ) ) THEN
         CALL PSOLVEJACTRANS( X, B )
      ELSE
         WRITE(*,*) 'IN PSOLVE: UNKNOWN PRECONDITIONER', CURPFORM, 
     $              ' QUITTING' 
         STOP
      ENDIF
*
      RETURN
*
      END
*     ================================================================
*
      SUBROUTINE PSOLVENONE( X, B )
*
*     .. Array Arguments ..
      DOUBLE PRECISION X( * ), B( * )
*     ..
*     .. Common Blocks ..
*
      INTEGER             N, LDA
      COMMON            / MATDIM / N, LDA
*
*  Purpose
*  =======
*
*  This PSOLVE is for the unpreconditioned version, i.e. just does
*  a vector copy ( B to X ) then returns.
*
*  Arguments
*  =========
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (output) DOUBLE PRECISION array, dimension N.
*          Set to solution on output.
*
*  BLAS:  DCOPY
*  ============================================================
*
*     .. Executable Statements ..
*
      CALL DCOPY( N, B, 1, X, 1 )

*
      RETURN
*
*     End of PSolveNone
*
      END
*
*     =====================================================
      SUBROUTINE PSOLVENONETRANS( X, B )
*
*     .. Array Arguments ..
*     ..
      DOUBLE PRECISION X( * ), B( * )
*     ..
*     .. Common Blocks ..
*
      INTEGER             N, LDA
      COMMON            / MATDIM / N, LDA
*     ..
*
*  Purpose
*  =======
*
*  This PSOLVE is for the unpreconditioned version, i.e. just does
*  a vector copy ( B to X ) then returns.
*
*  Arguments
*  =========
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (output) DOUBLE PRECISION array, dimension N.
*          Set to solution on output.
*
*  BLAS:  DCOPY
*  ============================================================
*
*     .. Executable Statements ..
*
      CALL DCOPY( N, B, 1, X, 1 )
*
      RETURN
*
*     End of PSolve
*
      END
*
*     ===========================================================
      SUBROUTINE PSOLVEJAC( X, B )
*
*     .. Array Delcarations ..
      DOUBLE PRECISION    X( * ), B( * )
*     ..
*     .. Common Blocks ..
*
      INTEGER             MAXDIM, MAXDIM2
      PARAMETER         ( MAXDIM = 200, MAXDIM2 = 40000 )
*
      INTEGER             N, LDA
      DOUBLE PRECISION    A, M
*
      COMMON            / SYSTEM / A( MAXDIM2 ), M( MAXDIM ),
     $                  / MATDIM / N, LDA
*
*  Purpose
*  =======
*
*  PSOLVE solves the linear system Mx = b where matrix M has 
*  is diagonal.
*
*  Arguments
*  =========
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (output) DOUBLE PRECISION array, dimension N.
*          Set to solution on output.
*  ============================================================
*
*     .. Local Scalars ..
*
      INTEGER          I
*
*     .. Executable Statements ..
*
      DO 10 I = 1, N
         X( I ) = B( I ) / M( I )
   10 CONTINUE
*
      RETURN
*
*     End of PSolveJac
*
      END
*
*     =========================================================
      SUBROUTINE PSOLVEJACTRANS( X, B )
*
*     .. Array Delcarations ..
*
      DOUBLE PRECISION X( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PSOLVETRANS solves the linear system Mx = b where matrix M has
*  is diagonal. Since this is the same as the non-transpose version,
*  this routine is actual just a mask to PSOLVEJAC.
*
*  Arguments
*  =========
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (output) DOUBLE PRECISION array, dimension N.
*          Set to solution on output.
*  ============================================================
*
*     .. Executable Statements ..
*
      CALL PSOLVEJAC( X, B )
*
      RETURN
*
*     End of PSolveJacTrans
*
      END
*     ================================================================
*     Following are the solvers for QMR, allowing left and right
*     preconditioning.
*     ================================================================
      SUBROUTINE PSOLVEQ( X, B, WHICH )
*
*     .. Argument Declarations ..
      DOUBLE PRECISION    X( * ), B( * )
      CHARACTER           WHICH*4
*     ..
*     .. Local Scalars ..
      LOGICAL             LSAME
*     ..
      CHARACTER           CURPFORM*5
      COMMON            / FORMS  / CURPFORM
*
*     .. Executable Statements ..
*
      IF ( LSAME( CURPFORM, 'IDENT' ) ) THEN
         CALL PSOLVENONE( X, B )
      ELSE IF ( LSAME( CURPFORM, 'JACBI' ) ) THEN
         IF( LSAME( WHICH, 'LEFT' ) ) THEN
            CALL PSOLVEJAC( X, B )
         ELSE
            CALL PSOLVENONE( X, B )
         ENDIF
      ELSE
         WRITE(*,*) 'IN PSOLVEQ: UNKNOWN PRECONDITIONER', CURPFORM, 
     $              ' QUITTING' 
         STOP
      ENDIF
*
      RETURN
*
      END
*     ================================================================
*
      SUBROUTINE PSOLVETRANSQ( X, B, WHICH )
*
*     .. Argument Declarations ..
      DOUBLE PRECISION X( * ), B( * )
      CHARACTER        WHICH*4
*     ..
*     .. Local Scalars ..
      LOGICAL          LSAME
*     ..
*     .. Common Blocks ..
      CHARACTER           CURPFORM*5
      COMMON            / FORMS  / CURPFORM
*
*     .. Executable Statements ..
*
      IF ( LSAME( CURPFORM, 'IDENT' ) ) THEN
         CALL PSOLVENONE( X, B )
      ELSE IF ( LSAME( CURPFORM, 'JACBI' ) ) THEN
         IF( LSAME( WHICH, 'LEFT' ) ) THEN
            CALL PSOLVEJAC( X, B )
         ELSE
            CALL PSOLVENONE( X, B )
         ENDIF
      ELSE
         WRITE(*,*) 'IN PSOLVEQ: UNKNOWN PRECONDITIONER', CURPFORM, 
     $              ' QUITTING' 
         STOP
      ENDIF
*
      RETURN
*
      END
