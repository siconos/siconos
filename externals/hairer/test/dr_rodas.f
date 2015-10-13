C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR ROSENBROCK CODE RODAS ON VAN DER POL
C * * * * * * * * * * * * * * * * * * * * * * * * *
c link dr_rodas rodas decsol dc_decsol
c link dr_rodas rodas lapack lapackc dc_lapack
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RODAS (FULL JACOBIAN)
        PARAMETER (ND=2,LWORK=2*ND*ND+14*ND+20,LIWORK=ND+20)
C --- DECLARATIONS
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
        EXTERNAL FVPOL,JVPOL,SOLOUT
C --- DIMENSION OF THE SYSTEM
        RPAR=1.D-6
        N=2
C --- PROBLEM IS AUTONOMOUS
        IFCN=0
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A FULL MATRIX
        MLJAC=N
C --- DIFFERENTIAL EQUATION IS IN EXPLICIT FORM
        IMAS=0
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=2.0D0
        Y(2)=-0.66D0
C --- ENDPOINT OF INTEGRATION
        XEND=2.0D0
C --- REQUIRED TOLERANCE
        RTOL=1.0D-4
        ATOL=1.0D-6*RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- SET DEFAULT VALUES
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.D0
C --- CALL OF THE SUBROUTINE RODAS
        CALL RODAS(N,FVPOL,IFCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JVPOL,IJAC,MLJAC,MUJAC,FVPOL,IDFX,
     &                  FVPOL,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
        WRITE (6,99) X,Y(1),Y(2)
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10)
C --- PRINT STATISTICS
        WRITE (6,90) RTOL
 90     FORMAT('       rtol=',D8.2)
        WRITE (6,91) (IWORK(J),J=14,20)
 91     FORMAT(' fcn=',I5,' jac=',I4,' step=',I4,
     &        ' accpt=',I4,' rejct=',I3,' dec=',I4,
     &        ' sol=',I5)
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
        COMMON /INTERN/XOUT
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),Y(2),NR-1
           XOUT=0.2D0
        ELSE
           IF (X.GE.XOUT) THEN
              Y1=CONTRO(1,XOUT,CONT,LRC)
              Y2=CONTRO(2,XOUT,CONT,LRC)
              WRITE (6,99) XOUT,Y1,Y2,NR-1
              XOUT=XOUT+0.2D0
           END IF
        END IF
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10,'    NSTEP =',I4)
        RETURN
        END
C
C
        SUBROUTINE FVPOL(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF VAN DER POL'S EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        F(1)=Y(2)
        F(2)=((1-Y(1)**2)*Y(2)-Y(1))/RPAR
        RETURN
        END
C
C
        SUBROUTINE JVPOL(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF VAN DER POL'S EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(LDFY,N)
        DFY(1,1)=0.0D0
        DFY(1,2)=1.0D0
        DFY(2,1)=(-2.0D0*Y(1)*Y(2)-1.0D0)/RPAR
        DFY(2,2)=(1.0D0-Y(1)**2)/RPAR
        RETURN
        END
