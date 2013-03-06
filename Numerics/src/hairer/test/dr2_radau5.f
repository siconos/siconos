C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR RADAU5 AT THE AMPLIFIER PROBLEM
C * * * * * * * * * * * * * * * * * * * * * * * * *
c link dr2_radau5 radau5 decsol dc_decsol
c link dr2_radau5 radau5 lapack lapackc dc_lapack
        IMPLICIT REAL*8 (A-H,O-Z)
C --- PARAMETERS FOR RADAU5 (BANDED JACOBIAN)
        PARAMETER (ND=8,LJAC=4,LMAS=3,LE=5)
        PARAMETER (LWORK=ND*(LJAC+LMAS+3*LE+12)+20,LIWORK=3*ND+20)
C --- DECLARATIONS
        DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),RPAR(15)
        EXTERNAL FAMPL,JBAMPL,BBAMPL,SOLOUT
C --- DATA FOR THE PROBLEM
        UE=0.1D0
          RPAR(1)=UE
        UB=6.0D0
          RPAR(2)=UB
        UF=0.026D0
          RPAR(3)=UF
        ALPHA=0.99D0
          RPAR(4)=ALPHA
        BETA=1.0D-6
          RPAR(5)=BETA
        R0=1000.0D0
          RPAR(6)=R0
        R1=9000.0D0
          RPAR(7)=R1
        R2=9000.0D0
          RPAR(8)=R2
        R3=9000.0D0
          RPAR(9)=R3
        R4=9000.0D0
          RPAR(10)=R4
        R5=9000.0D0
          RPAR(11)=R5
        R6=9000.0D0
          RPAR(12)=R6
        R7=9000.0D0
          RPAR(13)=R7
        R8=9000.0D0
          RPAR(14)=R8
        R9=9000.0D0
          RPAR(15)=R9
C --- DIMENSION OF THE SYSTEM
        N=8
C --- COMPUTE THE JACOBIAN ANALYTICALLY
        IJAC=1
C --- JACOBIAN IS A BANDED MATRIX (LOWER AND UPPER BAND-WIDTHS ARE 1, 2)
        MLJAC=1
        MUJAC=2
C --- DIFFERENTIAL EQUATION IS IN IMPLICIT FORM
        IMAS=1
        MLMAS=1
        MUMAS=1
C --- OUTPUT ROUTINE IS USED DURING INTEGRATION
        IOUT=1
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=0.D0
        Y(2)=UB-Y(1)*R8/R9
        Y(3)=UB/(R6/R5+1.D0)
        Y(4)=UB/(R6/R5+1.D0)
        Y(5)=UB
        Y(6)=UB/(R2/R1+1.D0)
        Y(7)=UB/(R2/R1+1.D0)
        Y(8)=0.D0
C --- ENDPOINT OF INTEGRATION
        XEND=0.05D0
C --- REQUIRED TOLERANCE
        RTOL=1.0D-5
        ATOL=1.0D-6*RTOL
        ITOL=0
C --- INITIAL STEP SIZE
        H=1.0D-6
C --- SET DEFAULT VALUES
        DO 10 I=1,20
        IWORK(I)=0
  10    WORK(I)=0.D0
C --- CALL OF THE SUBROUTINE RADAU5
        CALL RADAU5(N,FAMPL,X,Y,XEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JBAMPL,IJAC,MLJAC,MUJAC,
     &                  BBAMPL,IMAS,MLMAS,MUMAS,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
        WRITE (6,99) X,Y(1),Y(2)
 99     FORMAT(1X,'X =',F7.4,'    Y =',2E18.10)
C --- PRINT STATISTICS
        WRITE (6,90) RTOL
 90     FORMAT('       tol=',D8.2)
        WRITE (6,91) (IWORK(J),J=14,20)
 91     FORMAT(' fcn=',I5,' jac=',I4,' step=',I4,
     &        ' accpt=',I4,' rejct=',I3,' dec=',I4,
     &        ' sol=',I5)
        STOP
        END
C
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
C --- BY USING "CONTR5"
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
        COMMON /INTERN/XOUT
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),Y(2),NR-1
           XOUT=0.0025D0
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
C --- CONTINUOUS OUTPUT FOR RADAU5
              WRITE (6,99) XOUT,CONTR5(1,XOUT,CONT,LRC),
     &                     CONTR5(2,XOUT,CONT,LRC),NR-1
              XOUT=XOUT+0.0025D0
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F7.4,'    Y =',2E18.10,'    NSTEP =',I4)
        RETURN
        END
C
C
        SUBROUTINE FAMPL(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF THE AMPLIFIER PROBLEM
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 Y(N),F(N),RPAR(15)
        UE=RPAR(1)
        UB=RPAR(2)
        UF=RPAR(3)
        ALPHA=RPAR(4)
        BETA=RPAR(5)
        R0=RPAR(6)
        R1=RPAR(7)
        R2=RPAR(8)
        R3=RPAR(9)
        R4=RPAR(10)
        R5=RPAR(11)
        R6=RPAR(12)
        R7=RPAR(13)
        R8=RPAR(14)
        R9=RPAR(15)
        W=2.D0*3.141592654D0*100.D0
        UET=UE*DSIN(W*X)
        FAC1=BETA*(DEXP((Y(4)-Y(3))/UF)-1.D0)
        FAC2=BETA*(DEXP((Y(7)-Y(6))/UF)-1.D0)
        F(1)=Y(1)/R9
        F(2)=(Y(2)-UB)/R8+ALPHA*FAC1
        F(3)=Y(3)/R7-FAC1
        F(4)=Y(4)/R5+(Y(4)-UB)/R6+(1.D0-ALPHA)*FAC1
        F(5)=(Y(5)-UB)/R4+ALPHA*FAC2
        F(6)=Y(6)/R3-FAC2
        F(7)=Y(7)/R1+(Y(7)-UB)/R2+(1.D0-ALPHA)*FAC2
        F(8)=(Y(8)-UET)/R0
        RETURN
        END
C
C
        SUBROUTINE JBAMPL(N,X,Y,DFY,LDFY,RPAR,IPAR)
C --- JACOBIAN OF THE AMPLIFIER PROBLEM
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 Y(N),DFY(LDFY,N),RPAR(15)
        UE=RPAR(1)
        UB=RPAR(2)
        UF=RPAR(3)
        ALPHA=RPAR(4)
        BETA=RPAR(5)
        R0=RPAR(6)
        R1=RPAR(7)
        R2=RPAR(8)
        R3=RPAR(9)
        R4=RPAR(10)
        R5=RPAR(11)
        R6=RPAR(12)
        R7=RPAR(13)
        R8=RPAR(14)
        R9=RPAR(15)
        FAC14=BETA*DEXP((Y(4)-Y(3))/UF)/UF
        FAC27=BETA*DEXP((Y(7)-Y(6))/UF)/UF
        DO 1 J=1,8
        DFY(1,J)=0.D0
        DFY(2,J)=0.D0
   1    DFY(4,J)=0.D0
        DFY(3,1)=1.D0/R9
        DFY(3,2)=1.D0/R8
        DFY(2,3)=-ALPHA*FAC14
        DFY(1,4)=ALPHA*FAC14
        DFY(3,3)=1.D0/R7+FAC14
        DFY(2,4)=-FAC14
        DFY(3,4)=1.D0/R5+1.D0/R6+(1.D0-ALPHA)*FAC14
        DFY(4,3)=-(1.D0-ALPHA)*FAC14
        DFY(3,5)=1.D0/R4
        DFY(2,6)=-ALPHA*FAC27
        DFY(1,7)=ALPHA*FAC27
        DFY(3,6)=1.D0/R3+FAC27
        DFY(2,7)=-FAC27
        DFY(3,7)=1.D0/R1+1.D0/R2+(1.D0-ALPHA)*FAC27
        DFY(4,6)=-(1.D0-ALPHA)*FAC27
        DFY(3,8)=1.D0/R0
        RETURN
        END
C
        SUBROUTINE BBAMPL(N,B,LB,RPAR,IPAR)
C --- MATRIX "M" FOR THE AMPLIFIER PROBLEM
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 B(LB,N),RPAR(15)
        DO 1 I=1,8
        B(1,I)=0.D0
   1    B(3,I)=0.D0
        C1=1.D-6
        C2=2.D-6
        C3=3.D-6
        C4=4.D-6
        C5=5.D-6
C
        B(2,1)=-C5
        B(1,2)=C5
        B(3,1)=C5
        B(2,2)=-C5
        B(2,3)=-C4
        B(2,4)=-C3
        B(1,5)=C3
        B(3,4)=C3
        B(2,5)=-C3
        B(2,6)=-C2
        B(2,7)=-C1
        B(1,8)=C1
        B(3,7)=C1
        B(2,8)=-C1
        RETURN
        END
