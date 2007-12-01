*
*     This file contains routines for forming matrices that result from
*     a 5- or 7-point discretization of elliptic PDEs with Dirichlet
*     boundary conditions, and a consistent mass matrix "Wathen".
*
*        GEN57PT and GETSTEN are from SPARSEKIT. They actually form the
*        row compressed matrix.
*
*        COEFF provides the functions for computing the coefficients 
*        of the PDE.
*
*        Finally, for testing the iterative templates, COMP2DENSE converts
*        the row compressed matrix to dense form.
*
*     =================================================================
      SUBROUTINE GEN57PT( NX, NY, NZ, A, INDX, POINTR, AFUN, BFUN, CFUN,
     $                    DFUN, EFUN, FFUN, GFUN )
*
*     .. Scalar Arguments ..
      INTEGER            NX, NY, NZ
*     ..
*     .. Array Arguments ..
      INTEGER            INDX( * ), POINTR( * )
      DOUBLE PRECISION   A( * )
*     ..
*     .. External Functions ..
      EXTERNAL           AFUN, BFUN, CFUN, DFUN, EFUN, FFUN, GFUN
*
*  Purpose
*  =======
*
*  Adapted/altered from SPARSEKIT
*
*  This subroutine computes the sparse matrix in row compressed
*  format for the elliptic operator
*
*  L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) + 
*        delx ( d u ) + dely (e u) + delz( f u ) + g u 
*
*  with Dirichlet Boundary conditions, on a rectangular 1-D, 
*  2-D or 3-D grid using centered difference schemes.
* 
*  The functions a, b, ..., g are known through the
*  subroutines  afun, bfun, ..., gfun.
*  Note that to obtain the correct matrix, any function that is not
*  needed should be set to zero. For example for two-dimensional
*  problems, nz should be set to 1 and the functions cfun and ffun
*  should be zero functions. 
*
*  Uses natural ordering, first x direction, then y, then z
*  mesh size h is uniform and determined by grid points 
*  in the x-direction.
*
*  Arguments
*  =========
*
*  NX     (input) INTEGER
*         Number of points in X direction.
*
*  NY     (input) INTEGER
*         Number of points in Y direction.
*
*  NZ     (input) INTEGER
*         Number of points in Z direction.
*
*  A,     (output) DOUBLE PRECISION array.
*         Nonzero elements of the matrix. Stored in row compressed form.
* 
*  INDX   (output) INTEGER array.
*         Column index of matrix element.
* 
*  POINTR (output) INTEGER array.
*         Each element = P+1, where P is the number of nonzero elements
*         in the preceding rows of the matrix.
*
*  AFUN,
*  BFUN,
*  CFUN,
*  DFUN,
*  EFUN,
*  FFUN,
*  GFUN   (external subroutine) 
*         The user must supply the functions for computing the coefficients
*         of the PDE.
*
*  Description of the STENCIL:
*
*     stencil [1:7] has the following meaning:
*
*        center point = stencil(1)
*        west point   = stencil(2)
*        east point   = stencil(3)
*        south point  = stencil(4)
*        north point  = stencil(5)
*        front point  = stencil(6) 
*        back point   = stencil(7)
*
*                           st(5)
*                            |
*                            |  
*                            |
*                            |          .st(7)
*                            |     .
*                            | . 
*         st(2) ----------- st(1) ---------- st(3)
*                       .    |
*                   .        |
*               .            |
*            st(6)           |
*                            |
*                            |
*                           st(4)
*
*     ===============================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            KX, KY, KZ, IX, IY, IZ, IEDGE, NODE
      DOUBLE PRECISION   H, STENCIL( 7 )
*
*     .. Executable Statements ..
*
*     Initializations
*
      H     = ONE / (NX+1)
      KX    = 1
      KY    = NX
      KZ    = NX*NY
      IEDGE = 1
      NODE  = 1
*
      DO 30 IZ = 1, NZ
         DO 20 IY = 1, NY
            DO 10 IX = 1, NX
               POINTR(NODE) = IEDGE
*
*              Get stencil.
*
               CALL GETSTEN( NY, NZ, IX, IY, IZ, STENCIL, H,
     $                       AFUN, BFUN, CFUN, DFUN, EFUN, FFUN, GFUN )
*
*              West.
*
               IF ( IX.GT.1 ) THEN
                  INDX(IEDGE)=NODE-KX
                  A(IEDGE) = STENCIL(2)
                  IEDGE=IEDGE + 1
               ENDIF
*
*              South.
*
               IF ( IY.GT.1 ) THEN
                  INDX(IEDGE)=NODE-KY
                  A(IEDGE) = STENCIL(4)
                  IEDGE=IEDGE + 1
               END IF
*
*              Front Plane.
*
               IF ( IZ.GT.1 ) THEN
                  INDX(IEDGE)=NODE-KZ
                  A(IEDGE) = STENCIL(6)
                  IEDGE=IEDGE + 1
               ENDIF
*
*              Center node.
*
               INDX(IEDGE) = NODE
               A(IEDGE) = STENCIL(1)
               IEDGE = IEDGE + 1
*
*              Upper part.
*
*              East.
*
               IF ( IX.LT.NX ) THEN
                  INDX(IEDGE)=NODE+KX
                  A(IEDGE) = STENCIL(3)
                  IEDGE=IEDGE + 1
               END IF
               IF ( IY.LT.NY ) THEN
                  INDX(IEDGE)=NODE+KY
                  A(IEDGE) = STENCIL(5)
                  IEDGE=IEDGE + 1
               END IF
*
*              Back plane.
*
               IF ( IZ.LT.NZ ) THEN
                  INDX(IEDGE)=NODE+KZ
                  A(IEDGE) = STENCIL(7)
                  IEDGE=IEDGE + 1
               END IF
*
*              Next node.
*
               NODE=NODE+1
*
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
*
      POINTR(NODE)=IEDGE
*
      RETURN
*
*     -- End of GEN57PT
*
      END
*     ===============================================================
      SUBROUTINE GETSTEN( NY, NZ, KX, KY, KZ, STENCIL, H,
     $                    AFUN, BFUN, CFUN, DFUN, EFUN, FFUN, GFUN )
*
*     .. Argument Declarations ..
      INTEGER            NY, NZ, KX, KY, KZ
      DOUBLE PRECISION   H, AFUN, BFUN, CFUN, DFUN, EFUN, FFUN, GFUN,
     $                   STENCIL( * )
*     ..
*     .. External Functions ..
      EXTERNAL           AFUN, BFUN, CFUN, DFUN, EFUN, FFUN, GFUN
*
*  Purpose
*  =======
*
*  This subroutine calcultes the correct stencil values for
*  elliptic operator
*
*     L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) + 
*           delx ( d u ) + dely (e u) + delz( f u ) + g u.
*
*  For 2-D problems the discretization formula that is used is:
*      
*  h**2 * Lu == a(i+1/2,j)*{u(i+1,j) - u(i,j)} +
*               a(i-1/2,j)*{u(i-1,j) - u(i,j)} + 
*               b(i,j+1/2)*{u(i,j+1) - u(i,j)} +
*               b(i,j-1/2)*{u(i,j-1) - u(i,j)} + 
*              (h/2)*d(i,j)*{u(i+1,j) - u(i-1,j)} +
*              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} + 
*              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} + 
*              (h**2)*g(i,j)*u(i,j) 
*
*  ===================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF
      PARAMETER        ( ZERO = 0.0D+0, HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            K
      DOUBLE PRECISION   HHALF,CNTR, X, Y, Z, COEFF
*
*     .. Executable Statements ..
*
      DO 10 K=1,7
         STENCIL(K) = ZERO
  10  CONTINUE
*
      HHALF = H*HALF
      X = H*KX
      Y = H*KY
      Z = H*KZ
      CNTR = ZERO
*
*     Differentiation w.r.t. X.
*
      COEFF = AFUN( X+HHALF, Y, Z )
      STENCIL(3) = STENCIL(3) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = AFUN( X-HHALF,Y,Z )
      STENCIL(2) = STENCIL(2) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = DFUN( X, Y, Z )*HHALF
      STENCIL(3) = STENCIL(3) + COEFF
      STENCIL(2) = STENCIL(2) - COEFF
      IF (NY.LE.1) GOTO 99
*
*     Differentiation w.r.t. Y.
*
      COEFF = BFUN( X, Y+HHALF, Z )
      STENCIL(5) = STENCIL(5) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = BFUN( X, Y-HHALF, Z )
      STENCIL(4) = STENCIL(4) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = EFUN( X, Y, Z )*HHALF
      STENCIL(5) = STENCIL(5) + COEFF
      STENCIL(4) = STENCIL(4) - COEFF
      IF ( NZ.LE.1) GOTO 99
*
*     Differentiation w.r.t. Z.
*
      COEFF = CFUN( X, Y, Z+HHALF )
      STENCIL(7) = STENCIL(7) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = CFUN( X, Y, Z-HHALF )
      STENCIL(6) = STENCIL(6) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = FFUN( X, Y,Z )*HHALF
      STENCIL(7) = STENCIL(7) + COEFF
      STENCIL(6) = STENCIL(6) - COEFF
*
*     Discretization of product by G.
*
 99   COEFF = GFUN( X, Y, Z )
      STENCIL(1) = H*H*COEFF - CNTR
*
      RETURN
*
      END
*     =============================================================
*     Below are some functions for computing the value of the
*     coefficients.
*     =============================================================
      DOUBLE PRECISION FUNCTION ZEROFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      DOUBLE PRECISION   X, Y, Z
*
*     Purpose: Function to return ZERO.
*     =======
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*
*     .. Executable Statements ..
*
      ZEROFUN = ZERO
*
*     RETURN
*
      END
*
*     =============================================================
      DOUBLE PRECISION FUNCTION ONEFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      DOUBLE PRECISION   X, Y, Z
*
*     Purpose: Function to return ONE.
*     =======
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0D+0 )
*
*     .. Executable Statements ..
*
      ONEFUN = ONE
*
*     RETURN
*
      END
*
*     =============================================================
      DOUBLE PRECISION FUNCTION NEGONEFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      DOUBLE PRECISION   X, Y, Z
*
*     Purpose: Function to return -ONE.
*     =======
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE = 1.0D+0 )
*
*     .. Executable Statements ..
*
      NEGONEFUN = -ONE
*
*     RETURN
*
      END
*
*     =============================================================
      DOUBLE PRECISION FUNCTION THOUSFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      DOUBLE PRECISION   X, Y, Z
*
*     Purpose: Function to return 1000.
*     =======
*
*     .. Executable Statements ..
*
      THOUSFUN = 1000.0D+0
*
*     RETURN
*
      END
*
*     =============================================================
      DOUBLE PRECISION FUNCTION TEN5X2FUN( X, Y, Z )
*
*     .. Argument Declarations ..
      DOUBLE PRECISION   X, Y, Z
*
*     Purpose: Function to return 10 * X^2.
*     =======
*
*     .. Executable Statements ..
*
      TEN5X2FUN = 10.0D+5 * X * X
*
*     RETURN
*
      END
*
*     =============================================================
      DOUBLE PRECISION FUNCTION THOUSXFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      DOUBLE PRECISION   X, Y, Z
* 
*     Purpose: Evaluates the coefficient function
*     =======
*
*     .. Parameter ..
      DOUBLE PRECISION   NATX
      PARAMETER        ( NATX = 2.7182817459106D+0 )
*
*     .. Executable Statements ..
*
      THOUSXFUN = 1000.D+0 * ( NATX**( X*Y*Z ) )
*
*     RETURN
*
      END
*
*     =============================================================
      DOUBLE PRECISION FUNCTION NEGTHOUSXFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      DOUBLE PRECISION   X, Y, Z
*
*     Purpose: Evaluates the coefficient function
*     =======
*
*     .. Parameter ..
      DOUBLE PRECISION   NATX
      PARAMETER        ( NATX = 2.7182817459106D+0 )
*
*     .. Executable Statements ..
*
      NEGTHOUSXFUN = -1000.D+0 * ( NATX**( X*Y*Z ) )
*
*     RETURN
*
      END
*
*     =============================================================
      DOUBLE PRECISION FUNCTION HENKFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      DOUBLE PRECISION   X, Y, Z
*
*     Purpose: Evaluates the derivative of the above coefficient function
*     =======
*
*     .. Parameter ..
      DOUBLE PRECISION   NATX
      PARAMETER        ( NATX = 2.7182817459106D+0 )
*
*     .. Executable Statements ..
*
      HENKFUN = 20 * ( NATX**( 3.5D+0 * ( X**2 + Y**2 ) ) )
*
*     RETURN
*
      END
*
*     =============================================================
      DOUBLE PRECISION FUNCTION HENKDFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      DOUBLE PRECISION   X, Y, Z
*
*     Purpose: Evaluates the derivative of the above coefficient function.
*     =======
*
*     .. Parameter ..
      DOUBLE PRECISION   NATX
      PARAMETER        ( NATX = 2.7182817459106D+0 )
*
*     .. Executable Statements ..
*
      HENKDFUN = 70 * X * ( NATX**( 3.5D+0 * ( X**2 + Y**2 ) ) )
*
*     RETURN
*
      END
*
*     ===============================================================
      SUBROUTINE COMP2DENSE( ASPARSE, POINTR, INDX, N, ADENSE, LDA,
     $                       FLAG, INFO  )
*
      INTEGER            N, LDA, INFO, POINTR( * ), INDX( * )
      DOUBLE PRECISION   ASPARSE( * ), ADENSE( LDA,* )
      CHARACTER          FLAG*3
*
*     Convert sparse matrix storage to dense.
*
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
      INTEGER            I, J
      LOGICAL            LSAME, LSAMEN
*
*     .. Executable Statements ..
*
      INFO = 0
      IF ( N.LE.0 ) THEN
         INFO = -1
      ELSE IF ( N.GT.LDA ) THEN
         INFO = -2
      ELSE IF ( .NOT.LSAMEN( 3, FLAG,'ROW' ).OR.
     $          .NOT.LSAMEN( 3, FLAG,'ROW' ) ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
      DO 20 J = 1, N
         DO 10 I = 1, N
            ADENSE( I,J ) = ZERO
   10    CONTINUE
   20 CONTINUE
*
      IF ( LSAME( FLAG,'ROW' ) ) THEN
         DO 40 I = 1, N
            DO 30 J = POINTR( I ), POINTR( I+1 ) - 1
               ADENSE( I,INDX( J ) ) = ASPARSE( J )
   30       CONTINUE
   40    CONTINUE
      ELSE IF ( LSAME( FLAG,'COL' ) ) THEN
         DO 60 J = 1, N
            DO 50 I = POINTR( J ), POINTR( J+1 ) - 1
               ADENSE( INDX( I ),J ) = ASPARSE( I ) 
   50       CONTINUE
   60    CONTINUE
      ENDIF
*
      RETURN
*
      END
*
*     ================================================================
      SUBROUTINE WATHEN( NX, NY, KK, N, A, LDA, WORK, LDW, INFO )
*
      INTEGER            NX, NY, KK, N, LDA, LDW, INFO
      DOUBLE PRECISION   A( LDA,* ), WORK( LDW,* )
*
*     Translated from the matlab version found on netlib.
*
*     A is a random N-by-N finite element matrix where 
*     N = 3*NX*NY + 2*NX + 2*NY + 1. A is precisely the "consistent 
*     mass matrix" for a regular NX-by-NY grid of 8-node (serendipity) 
*     elements in 2 space dimensions. A is symmetric positive definite 
*     for any (positive) values of the "density", RHO(NX,NY), which is 
*     chosen randomly in this routine. In particular, if D=DIAG(DIAG(A)), 
*     then 0.25 <= EIG(INV(D)*A) <= 4.5 for any positive integers NX and NY 
*     and any densities RHO(NX,NY). This diagonally scaled matrix is 
*     returned by WATHEN(NX,NY,1).
*
*     Reference: A.J.Wathen, DOUBLE PRECISIONistic eigenvalue bounds for 
*     the Galerkin mass matrix, IMA J. Numer. Anal., 7 (1987), pp. 449-457.
*
*     BEWARE - this is a sparse matrix, stored in -dense- form, and 
*              it quickly gets large!
*
*     .. Local Scalars ..
*
      INTEGER            I, J, E, EM, RHO, KROW, KCOL,
     $                   ISEED( 4 ), NN( 8 )
      DOUBLE PRECISION   RHOIT, DLARAN, ZERO, ONE
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0  )
*
*     .. Executable Statements ..
*
      INFO = 0
      N = 3*NX*NY + 2*NX + 2*NY + 1
      IF ( N.GT.LDA ) THEN
         WRITE(*,*) 'NOT ENOUGH ROOM ALLOCATED FOR WATHEN MATRIX'
         INFO = -1
         RETURN
      ELSE IF ( NX.LT.1 ) THEN
         INFO = -2
      ELSE IF ( NY.LT.1 ) THEN
         INFO = -3
      ELSE IF ( MAX( NX,NY ).GT.LDW ) THEN
         INFO = -4
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Alias workspace columns.
*
      E   = 1
      EM  = E + 8
      RHO = EM + 8
*
      CALL SET_E( WORK( 1,E ), LDW )
*
      DO 20 J = 1, N
         DO 10 I = 1, N
            A( I,J ) = ZERO
   10    CONTINUE
   20 CONTINUE
*
      ISEED( 1 ) = 304
      ISEED( 2 ) = 152
      ISEED( 3 ) = 2042
      ISEED( 4 ) = 77
      DO 40 J = 1, NY
         DO 30 I = 1, NX
            WORK( I,RHO+J-1 ) = 100*DLARAN( ISEED )
   30    CONTINUE
   40 CONTINUE
*
      DO 100 J = 1, NY
         DO 90 I = 1, NX
*
            NN(1) = 3*J*NX + 2*I + 2*J + 1
            NN(2) = NN(1)-1
            NN(3) = NN(2)-1
            NN(4) = (3*J-1)*NX+2*J+I-1
            NN(5) = 3*(J-1)*NX+2*I+2*J-3
            NN(6) = NN(5)+1
            NN(7) = NN(6)+1
            NN(8) = NN(4)+1
*
            RHOIT = WORK( I,RHO+J-1 )
            DO 60 KROW = 1, 8
               DO 50 KCOL = 1, 8
                  WORK( KROW,EM+KCOL-1 ) = RHOIT * WORK( KROW,E+KCOL-1 )
   50          CONTINUE
   60       CONTINUE
*
            DO 80 KROW = 1, 8
               DO 70 KCOL = 1, 8
                  A(NN(KROW),NN(KCOL)) = A(NN(KROW),NN(KCOL)) +
     $                                      WORK( KROW,EM+KCOL-1)
   70          CONTINUE
   80       CONTINUE 
*
   90    CONTINUE
  100 CONTINUE
*  
      IF ( KK.EQ.1 ) THEN
*
*        A = diag(diag(A)) \ A (the result being unit diagonal);
*
         DO 110 I = 1, J
            A( I,I ) = ONE
  110    CONTINUE 
      ENDIF
*
      RETURN
*
      END
*     ========================================================
      SUBROUTINE SET_E( E, LDE )
*
      INTEGER            I, J, LDE
      DOUBLE PRECISION   SCALE, E( LDE,* )
*
      E( 1,1 ) =  6.0D+0
      E( 2,1 ) = -6.0D+0
      E( 3,1 ) =  2.0D+0
      E( 4,1 ) = -8.0D+0
      E( 1,2 ) = -6.0D+0
      E( 2,2 ) = 32.0D+0
      E( 3,2 ) = -6.0D+0
      E( 4,2 ) = 20.0D+0
      E( 1,3 ) =  2.0D+0
      E( 2,3 ) = -6.0D+0
      E( 3,3 ) =  6.0D+0
      E( 4,3 ) = -6.0D+0
      E( 1,4 ) = -8.0D+0
      E( 2,4 ) = 20.0D+0
      E( 3,4 ) = -6.0D+0
      E( 4,4 ) = 32.0D+0
*
      E( 1,5 ) =  3.0D+0
      E( 2,5 ) = -8.0D+0
      E( 3,5 ) =  2.0D+0
      E( 4,5 ) = -6.0D+0
      E( 1,6 ) = -8.0D+0
      E( 2,6 ) = 16.0D+0
      E( 3,6 ) = -8.0D+0
      E( 4,6 ) = 20.0D+0
      E( 1,7 ) =  2.0D+0
      E( 2,7 ) = -8.0D+0
      E( 3,7 ) =  3.0D+0
      E( 4,7 ) = -8.0D+0
      E( 1,8 ) = -6.0D+0
      E( 2,8 ) = 20.0D+0
      E( 3,8 ) = -8.0D+0
      E( 4,8 ) = 16.0D+0
*
      DO 20 J = 1, 4
         DO 10 I = 5, 8 
            E( I,J ) = E( J,I ) 
   10    CONTINUE
   20 CONTINUE
*
      DO 40 J = 5, 8
         DO 30 I = 5, 8 
            E( I,J ) = E( I-4,J-4 )
   30    CONTINUE
   40 CONTINUE
*
      SCALE = 1.0D+0 / 45.0D+0
      DO 50 I = 1, 8
         CALL DSCAL( 8, SCALE, E( 1,I ), 1 )
   50 CONTINUE
*
      RETURN
*
      END
