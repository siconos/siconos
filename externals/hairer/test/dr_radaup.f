C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAUP AT VAN DER POL'S EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
c link dr_radaup radaup decsol dc_decsol
c link dr_radaup radaup lapack lapackc dc_lapack
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (FULL JACOBIAN)
        PARAMETER (ND=2,NS=7,LWORK=(NS+1)*ND*ND+(3*NS+3)*ND+20,
     &             LIWORK=(2+(NS-1)/2)*ND+20)
C --- DECLARATIONS
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK)
        EXTERNAL FVPOL,JVPOL,SOLOUT
        RPAR=1.0D-6
        DO IORD=3,7,2
C --- DIMENSION OF THE SYSTEM
        N=2
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
        ATOL=1.0D0*RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- SET DEFAULT VALUES
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.D0
        IWORK(11)=IORD
C --- CALL OF THE SUBROUTINE RADAU5 (OR SDIRK4)
        CALL RADAUP(N,FVPOL,X,Y,XEND,H,
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
        END DO
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
C --- BY USING "CONTRP"
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
        COMMON /INTERN/XOUT
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),Y(2),NR-1
           XOUT=0.2D0
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
C --- CONTINUOUS OUTPUT FOR RADAUP
              WRITE (6,99) XOUT,CONTRP(1,XOUT,CONT,LRC),
     &                     CONTRP(2,XOUT,CONT,LRC),NR-1
              XOUT=XOUT+0.2D0
              GOTO 10
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
