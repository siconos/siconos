C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR EXTRAPOLATION CODES SEULEX AT VAN DER POL
C * * * * * * * * * * * * * * * * * * * * * * * * *
clink dr_seulex seulex dc_decsol decsol
clink dr_seulex seulex dc_lapack lapack lapackc
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR SEULEX (FULL JACOBIAN)
        PARAMETER (ND=2,KM=12,KM2=2+KM*(KM+3)/2,NRDENS=ND)
        PARAMETER (LWORK=2*ND*ND+(KM+8)*ND+4*KM+20+KM2*NRDENS)
        PARAMETER (LIWORK=2*ND+KM+20+NRDENS)
C --- DECLARATIONS
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
        EXTERNAL FVPOL,JVPOL,SOLOUT
        RPAR=1.0D-6
C --- DIMENSION OF THE SYSTEM
        N=2
C --- PROBLEM IS AUTONOMOUS
        IFCN=1
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=0
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
        ATOL=RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- SET DEFAULT VALUES
        DO 10 I=1,20
  10    WORK(I)=0.D0
        DO 12 I=1,20
  12    IWORK(I)=0
        IWORK(6)=N
        DO I=1,NRDENS
          IWORK(20+I)=I
        END DO
C --- CALL OF THE SUBROUTINE SEULEX
        CALL SEULEX(N,FVPOL,IFCN,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JVPOL,IJAC,MLJAC,MUJAC,
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
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,RC,LRC,IC,LIC,N,
     &                     RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),RC(LRC),IC(LIC)
        COMMON /INTERN/XOUT
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),Y(2),NR-1
           XOUT=0.2D0
        ELSE
           IF (X.GE.XOUT) THEN
              WRITE (6,99) X,Y(1),Y(2),NR-1
c             write (6,*) CONTEX(1,XOLD,RC,LRC,IC,LIC)
              XOUT=MAX(XOUT+0.2D0,X)
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
