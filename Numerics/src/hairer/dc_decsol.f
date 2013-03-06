C ******************************************
C     VERSION OF SEPTEMBER 18, 1995
C ******************************************
C
      SUBROUTINE DECOMR(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &            M1,M2,NM1,FAC1,E1,LDE1,IP1,IER,IJOB,CALHES,IPHES)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),
     &          IP1(NM1),IPHES(N)
      LOGICAL CALHES
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO  I=1,N
            E1(I,J)=-FJAC(I,J)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=-FJAC(I,JM1)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
 45   MM=M1/M2
      DO J=1,M2
         DO I=1,NM1
            SUM=0.D0
            DO K=0,MM-1
               SUM=(SUM+FJAC(I,J+K*M2))/FAC1
            END DO
            E1(I,J)=E1(I,J)-SUM
         END DO
      END DO
      CALL DEC (NM1,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,J)
         END DO
         E1(MDIAG,J)=E1(MDIAG,J)+FAC1
      END DO
      CALL DECB (N,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,JM1)
         END DO
         E1(MDIAG,J)=E1(MDIAG,J)+FAC1
      END DO
  46  MM=M1/M2
      DO J=1,M2
         DO I=1,MBJAC
            SUM=0.D0
            DO K=0,MM-1
               SUM=(SUM+FJAC(I,J+K*M2))/FAC1
            END DO
            E1(I+MLE,J)=E1(I+MLE,J)-SUM
         END DO
      END DO
      CALL DECB (NM1,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E1(I,J)=-FJAC(I,J)
         END DO
         DO I=MAX(1,J-MUMAS),MIN(N,J+MLMAS)
            E1(I,J)=E1(I,J)+FAC1*FMAS(I-J+MBDIAG,J)
         END DO
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=-FJAC(I,JM1)
         END DO
         DO I=MAX(1,J-MUMAS),MIN(NM1,J+MLMAS)
            E1(I,J)=E1(I,J)+FAC1*FMAS(I-J+MBDIAG,J)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,J)
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            E1(IB,J)=E1(IB,J)+FAC1*FMAS(I,J)
         END DO
      END DO
      CALL DECB (N,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  14  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,JM1)
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            E1(IB,J)=E1(IB,J)+FAC1*FMAS(I,J)
         END DO
      END DO
      GOTO 46
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E1(I,J)=FMAS(I,J)*FAC1-FJAC(I,J)
         END DO
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=FMAS(I,J)*FAC1-FJAC(I,JM1)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      IF (CALHES) CALL ELMHES (LDJAC,N,1,N,FJAC,IPHES)
      CALHES=.FALSE.
      DO J=1,N-1
         J1=J+1
         E1(J1,J)=-FJAC(J1,J)
      END DO
      DO J=1,N
         DO I=1,J
            E1(I,J)=-FJAC(I,J)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
      CALL DECH(N,LDE1,E1,1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE DECOMR
C
C ***********************************************************
C
      SUBROUTINE DECOMC(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &            M1,M2,NM1,ALPHN,BETAN,E2R,E2I,LDE1,IP2,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),
     &          E2R(LDE1,NM1),E2I(LDE1,NM1),IP2(NM1)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E2R(I,J)=-FJAC(I,J)
            E2I(I,J)=0.D0
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
      CALL DECC (N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=-FJAC(I,JM1)
            E2I(I,J)=0.D0
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
  45  MM=M1/M2
      ABNO=ALPHN**2+BETAN**2
      ALP=ALPHN/ABNO
      BET=BETAN/ABNO
      DO J=1,M2
         DO I=1,NM1
            SUMR=0.D0
            SUMI=0.D0
            DO K=0,MM-1
               SUMS=SUMR+FJAC(I,J+K*M2)
               SUMR=SUMS*ALP+SUMI*BET
               SUMI=SUMI*ALP-SUMS*BET
            END DO
            E2R(I,J)=E2R(I,J)-SUMR
            E2I(I,J)=E2I(I,J)-SUMI
         END DO
      END DO
      CALL DECC (NM1,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            IMLE=I+MLE
            E2R(IMLE,J)=-FJAC(I,J)
            E2I(IMLE,J)=0.D0
         END DO
         E2R(MDIAG,J)=E2R(MDIAG,J)+ALPHN
         E2I(MDIAG,J)=BETAN
      END DO
      CALL DECBC (N,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E2R(I+MLE,J)=-FJAC(I,JM1)
            E2I(I+MLE,J)=0.D0
         END DO
         E2R(MDIAG,J)=E2R(MDIAG,J)+ALPHN
         E2I(MDIAG,J)=E2I(MDIAG,J)+BETAN
      END DO
  46  MM=M1/M2
      ABNO=ALPHN**2+BETAN**2
      ALP=ALPHN/ABNO
      BET=BETAN/ABNO
      DO J=1,M2
         DO I=1,MBJAC
            SUMR=0.D0
            SUMI=0.D0
            DO K=0,MM-1
               SUMS=SUMR+FJAC(I,J+K*M2)
               SUMR=SUMS*ALP+SUMI*BET
               SUMI=SUMI*ALP-SUMS*BET
            END DO
            IMLE=I+MLE
            E2R(IMLE,J)=E2R(IMLE,J)-SUMR
            E2I(IMLE,J)=E2I(IMLE,J)-SUMI
         END DO
      END DO
      CALL DECBC (NM1,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO  J=1,N
         DO  I=1,N
            E2R(I,J)=-FJAC(I,J)
            E2I(I,J)=0.D0
         END DO
      END DO
      DO J=1,N
         DO I=MAX(1,J-MUMAS),MIN(N,J+MLMAS)
            BB=FMAS(I-J+MBDIAG,J)
            E2R(I,J)=E2R(I,J)+ALPHN*BB
            E2I(I,J)=BETAN*BB
         END DO
      END DO
      CALL DECC(N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=-FJAC(I,JM1)
            E2I(I,J)=0.D0
         END DO
         DO I=MAX(1,J-MUMAS),MIN(NM1,J+MLMAS)
            FFMA=FMAS(I-J+MBDIAG,J)
            E2R(I,J)=E2R(I,J)+ALPHN*FFMA
            E2I(I,J)=E2I(I,J)+BETAN*FFMA
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            IMLE=I+MLE
            E2R(IMLE,J)=-FJAC(I,J)
            E2I(IMLE,J)=0.D0
         END DO
         DO I=MAX(1,MUMAS+2-J),MIN(MBB,MUMAS+1-J+N)
            IB=I+MDIFF
            BB=FMAS(I,J)
            E2R(IB,J)=E2R(IB,J)+ALPHN*BB
            E2I(IB,J)=BETAN*BB
         END DO
      END DO
      CALL DECBC (N,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  14  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E2R(I+MLE,J)=-FJAC(I,JM1)
            E2I(I+MLE,J)=0.D0
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            FFMA=FMAS(I,J)
            E2R(IB,J)=E2R(IB,J)+ALPHN*FFMA
            E2I(IB,J)=E2I(IB,J)+BETAN*FFMA
         END DO
      END DO
      GOTO 46
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            BB=FMAS(I,J)
            E2R(I,J)=BB*ALPHN-FJAC(I,J)
            E2I(I,J)=BB*BETAN
         END DO
      END DO
      CALL DECC(N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=ALPHN*FMAS(I,J)-FJAC(I,JM1)
            E2I(I,J)=BETAN*FMAS(I,J)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO J=1,N-1
         J1=J+1
         E2R(J1,J)=-FJAC(J1,J)
         E2I(J1,J)=0.D0
      END DO
      DO J=1,N
         DO I=1,J
            E2I(I,J)=0.D0
            E2R(I,J)=-FJAC(I,J)
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
      CALL DECHC(N,LDE1,E2R,E2I,1,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE DECOMC
C
C ***********************************************************
C
      SUBROUTINE SLVRAR(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,E1,LDE1,Z1,F1,IP1,IPHES,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),
     &          IP1(NM1),IPHES(N),Z1(N),F1(N)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
 48   CONTINUE
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            DO I=1,NM1
               IM1=I+M1
               Z1(IM1)=Z1(IM1)+FJAC(I,JKM)*SUM1
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,Z1(M1+1),IP1)
 49   CONTINUE
      DO I=M1,1,-1
         Z1(I)=(Z1(I)+Z1(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
  45  CONTINUE
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               Z1(IM1)=Z1(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM1
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,Z1(M1+1),IP1)
      GOTO 49
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            S1=S1-FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            S1=S1-FMAS(I-J+MBDIAG,J)*F1(J+M1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
      END DO
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S1=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            S1=S1-FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         DO J=1,N
            S1=S1-FMAS(I,J)*F1(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         DO J=1,NM1
            S1=S1-FMAS(I,J)*F1(J+M1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
      END DO
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         Z1(I)=Z1(I)-F1(I)*FAC1
      END DO
      DO MM=N-2,1,-1
          MP=N-MM
          MP1=MP-1
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 746
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
 746      CONTINUE
          DO I=MP+1,N
             Z1(I)=Z1(I)-FJAC(I,MP1)*Z1(MP)
          END DO
       END DO
       CALL SOLH(N,LDE1,E1,1,Z1,IP1)
       DO MM=1,N-2
          MP=N-MM
          MP1=MP-1
          DO I=MP+1,N
             Z1(I)=Z1(I)+FJAC(I,MP1)*Z1(MP)
          END DO
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 750
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
 750      CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVRAR
C
C ***********************************************************
C
      SUBROUTINE SLVRAI(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,ALPHN,BETAN,E2R,E2I,LDE1,Z2,Z3,
     &          F2,F3,CONT,IP2,IPHES,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),
     &          IP2(NM1),IPHES(N),Z2(N),Z3(N),F2(N),F3(N)
      DIMENSION E2R(LDE1,NM1),E2I(LDE1,NM1)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLC (N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
 48   ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=1,NM1
               IM1=I+M1
               Z2(IM1)=Z2(IM1)+FJAC(I,JKM)*SUM2
               Z3(IM1)=Z3(IM1)+FJAC(I,JKM)*SUM3
            END DO
         END DO
      END DO
      CALL SOLC (NM1,LDE1,E2R,E2I,Z2(M1+1),Z3(M1+1),IP2)
 49   CONTINUE
      DO I=M1,1,-1
         MPI=M2+I
         Z2I=Z2(I)+Z2(MPI)
         Z3I=Z3(I)+Z3(MPI)
         Z3(I)=(Z3I*ALPHN-Z2I*BETAN)/ABNO
         Z2(I)=(Z2I*ALPHN+Z3I*BETAN)/ABNO
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLBC (N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
  45  ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               IIMU=I+MUJAC+1-J
               Z2(IM1)=Z2(IM1)+FJAC(IIMU,JKM)*SUM2
               Z3(IM1)=Z3(IM1)+FJAC(IIMU,JKM)*SUM3
            END DO
         END DO
      END DO
      CALL SOLBC (NM1,LDE1,E2R,E2I,MLE,MUE,Z2(M1+1),Z3(M1+1),IP2)
      GOTO 49
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            JM1=J+M1
            BB=FMAS(I-J+MBDIAG,J)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLBC(N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=0.0D0
         S3=0.0D0
         DO J=1,N
            BB=FMAS(I,J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S2=0.0D0
         S3=0.0D0
         DO J=1,NM1
            JM1=J+M1
            BB=FMAS(I,J)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO MM=N-2,1,-1
          MP=N-MM
          MP1=MP-1
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 746
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 746      CONTINUE
          DO I=MP+1,N
             E1IMP=FJAC(I,MP1)
             Z2(I)=Z2(I)-E1IMP*Z2(MP)
             Z3(I)=Z3(I)-E1IMP*Z3(MP)
          END DO
       END DO
       CALL SOLHC(N,LDE1,E2R,E2I,1,Z2,Z3,IP2)
       DO MM=1,N-2
          MP=N-MM
          MP1=MP-1
          DO I=MP+1,N
             E1IMP=FJAC(I,MP1)
             Z2(I)=Z2(I)+E1IMP*Z2(MP)
             Z3(I)=Z3(I)+E1IMP*Z3(MP)
          END DO
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 750
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 750      CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVRAI
C
C ***********************************************************
C
      SUBROUTINE SLVRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,ALPHN,BETAN,E1,E2R,E2I,LDE1,Z1,Z2,Z3,
     &          F1,F2,F3,CONT,IP1,IP2,IPHES,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),
     &          E2R(LDE1,NM1),E2I(LDE1,NM1),IP1(NM1),IP2(NM1),
     &          IPHES(N),Z1(N),Z2(N),Z3(N),F1(N),F2(N),F3(N)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC (N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
 48   ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=1,NM1
               IM1=I+M1
               Z1(IM1)=Z1(IM1)+FJAC(I,JKM)*SUM1
               Z2(IM1)=Z2(IM1)+FJAC(I,JKM)*SUM2
               Z3(IM1)=Z3(IM1)+FJAC(I,JKM)*SUM3
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,Z1(M1+1),IP1)
      CALL SOLC (NM1,LDE1,E2R,E2I,Z2(M1+1),Z3(M1+1),IP2)
 49   CONTINUE
      DO I=M1,1,-1
         MPI=M2+I
         Z1(I)=(Z1(I)+Z1(MPI))/FAC1
         Z2I=Z2(I)+Z2(MPI)
         Z3I=Z3(I)+Z3(MPI)
         Z3(I)=(Z3I*ALPHN-Z2I*BETAN)/ABNO
         Z2(I)=(Z2I*ALPHN+Z3I*BETAN)/ABNO
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      CALL SOLBC (N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
  45  ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               FFJA=FJAC(I+MUJAC+1-J,JKM)
               Z1(IM1)=Z1(IM1)+FFJA*SUM1
               Z2(IM1)=Z2(IM1)+FFJA*SUM2
               Z3(IM1)=Z3(IM1)+FFJA*SUM3
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,Z1(M1+1),IP1)
      CALL SOLBC (NM1,LDE1,E2R,E2I,MLE,MUE,Z2(M1+1),Z3(M1+1),IP2)
      GOTO 49
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         J1B=MAX(1,I-MLMAS)
         J2B=MIN(NM1,I+MUMAS)
         DO J=J1B,J2B
            JM1=J+M1
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(JM1)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      CALL SOLBC(N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=1,N
            BB=FMAS(I,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=1,NM1
            JM1=J+M1
            BB=FMAS(I,J)
            S1=S1-BB*F1(JM1)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO MM=N-2,1,-1
          MP=N-MM
          MP1=MP-1
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 746
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 746      CONTINUE
          DO I=MP+1,N
             E1IMP=FJAC(I,MP1)
             Z1(I)=Z1(I)-E1IMP*Z1(MP)
             Z2(I)=Z2(I)-E1IMP*Z2(MP)
             Z3(I)=Z3(I)-E1IMP*Z3(MP)
          END DO
       END DO
       CALL SOLH(N,LDE1,E1,1,Z1,IP1)
       CALL SOLHC(N,LDE1,E2R,E2I,1,Z2,Z3,IP2)
       DO MM=1,N-2
          MP=N-MM
          MP1=MP-1
          DO I=MP+1,N
             E1IMP=FJAC(I,MP1)
             Z1(I)=Z1(I)+E1IMP*Z1(MP)
             Z2(I)=Z2(I)+E1IMP*Z2(MP)
             Z3(I)=Z3(I)+E1IMP*Z3(MP)
          END DO
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 750
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 750      CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVRAD
C
C ***********************************************************
C
      SUBROUTINE ESTRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          H,DD1,DD2,DD3,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,
     &          E1,LDE1,Z1,Z2,Z3,CONT,F1,F2,IP1,IPHES,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),IP1(NM1),
     &     SCAL(N),IPHES(N),Z1(N),Z2(N),Z3(N),F1(N),F2(N),Y0(N),Y(N)
      DIMENSION CONT(N),RPAR(1),IPAR(1)
      LOGICAL FIRST,REJECT
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
      HEE1=DD1/H
      HEE2=DD2/H
      HEE3=DD3/H
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
   1  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO  I=1,N
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  11  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
  48  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=1,NM1
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1)
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   2  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  12  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
  45  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   3  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  13  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   4  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  14  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 45
C
   5  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=1,N
            SUM=SUM+FMAS(I,J)*F1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  15  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=1,NM1
            SUM=SUM+FMAS(I,J)*F1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   6  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ------  THIS OPTION IS NOT PROVIDED
      RETURN
C
   7  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
      END DO
      DO MM=N-2,1,-1
         MP=N-MM
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 310
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 310     CONTINUE
         DO I=MP+1,N
            CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
         END DO
      END DO
      CALL SOLH(N,LDE1,E1,1,CONT,IP1)
      DO MM=1,N-2
         MP=N-MM
         DO I=MP+1,N
            CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
         END DO
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 440
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 440     CONTINUE
      END DO
C
C --------------------------------------
C
  77  CONTINUE
      ERR=0.D0
      DO  I=1,N
         ERR=ERR+(CONT(I)/SCAL(I))**2
      END DO
      ERR=MAX(SQRT(ERR/N),1.D-10)
C
      IF (ERR.LT.1.D0) RETURN
      IF (FIRST.OR.REJECT) THEN
          DO I=1,N
             CONT(I)=Y(I)+CONT(I)
          END DO
          CALL FCN(N,X,CONT,F1,RPAR,IPAR)
          NFCN=NFCN+1
          DO I=1,N
             CONT(I)=F1(I)+F2(I)
          END DO
          GOTO (31,32,31,32,31,32,33,55,55,55,41,42,41,42,41), IJOB
C ------ FULL MATRIX OPTION
  31      CONTINUE
          CALL SOL(N,LDE1,E1,CONT,IP1)
          GOTO 88
C ------ FULL MATRIX OPTION, SECOND ORDER
 41      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=1,NM1
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOL(NM1,LDE1,E1,CONT(M1+1),IP1)
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
         GOTO 88
C ------ BANDED MATRIX OPTION
 32      CONTINUE
         CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
         GOTO 88
C ------ BANDED MATRIX OPTION, SECOND ORDER
 42      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
          GOTO 88
C ------ HESSENBERG MATRIX OPTION
  33      CONTINUE
          DO MM=N-2,1,-1
             MP=N-MM
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 510
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 510         CONTINUE
             DO I=MP+1,N
                CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
             END DO
          END DO
          CALL SOLH(N,LDE1,E1,1,CONT,IP1)
          DO MM=1,N-2
             MP=N-MM
             DO I=MP+1,N
                CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
             END DO
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 640
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 640         CONTINUE
          END DO
C -----------------------------------
   88     CONTINUE
          ERR=0.D0
          DO I=1,N
             ERR=ERR+(CONT(I)/SCAL(I))**2
          END DO
          ERR=MAX(SQRT(ERR/N),1.D-10)
       END IF
       RETURN
C -----------------------------------------------------------
  55   CONTINUE
       RETURN
       END
C
C     END OF SUBROUTINE ESTRAD
C
C ***********************************************************
C
      SUBROUTINE ESTRAV(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          H,DD,FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,NS,NNS,
     &          E1,LDE1,ZZ,CONT,FF,IP1,IPHES,SCAL,ERR,
     &          FIRST,REJECT,FAC1,RPAR,IPAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E1(LDE1,NM1),IP1(NM1),
     &     SCAL(N),IPHES(N),ZZ(NNS),FF(NNS),Y0(N),Y(N)
      DIMENSION DD(NS),CONT(N),RPAR(1),IPAR(1)
      LOGICAL FIRST,REJECT
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
   1  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  11  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
  48  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=1,NM1
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1)
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   2  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  12  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
  45  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   3  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J)
         END DO
         FF(I+N)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  13  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO  I=1,M1
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO I=M1+1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J+M1)
         END DO
         IM1=I+M1
         FF(IM1+N)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   4  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J)
         END DO
         FF(I+N)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  14  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO  I=1,M1
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO I=M1+1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*FF(J+M1)
         END DO
         IM1=I+M1
         FF(IM1+N)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 45
C
   5  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,N
         SUM=0.D0
         DO J=1,N
            SUM=SUM+FMAS(I,J)*FF(J)
         END DO
         FF(I+N)=SUM
         CONT(I)=SUM+Y0(I)
      END DO
      CALL SOL (N,LDE1,E1,CONT,IP1)
      GOTO 77
C
  15  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO  I=1,M1
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO I=M1+1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I)=SUM/H
      END DO
      DO I=1,NM1
         SUM=0.D0
         DO J=1,NM1
            SUM=SUM+FMAS(I,J)*FF(J+M1)
         END DO
         IM1=I+M1
         FF(IM1+N)=SUM
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      GOTO 48
C
   6  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ------  THIS OPTION IS NOT PROVIDED
      RETURN
C
   7  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO  I=1,N
         SUM=0.D0
         DO K=1,NS
            SUM=SUM+DD(K)*ZZ(I+(K-1)*N)
         END DO
         FF(I+N)=SUM/H
         CONT(I)=FF(I+N)+Y0(I)
      END DO
      DO MM=N-2,1,-1
         MP=N-MM
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 310
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 310     CONTINUE
         DO I=MP+1,N
            CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
         END DO
      END DO
      CALL SOLH(N,LDE1,E1,1,CONT,IP1)
      DO MM=1,N-2
         MP=N-MM
         DO I=MP+1,N
            CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
         END DO
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 440
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
 440     CONTINUE
      END DO
C
C --------------------------------------
C
  77  CONTINUE
      ERR=0.D0
      DO  I=1,N
         ERR=ERR+(CONT(I)/SCAL(I))**2
      END DO
      ERR=MAX(SQRT(ERR/N),1.D-10)
C
      IF (ERR.LT.1.D0) RETURN
      IF (FIRST.OR.REJECT) THEN
          DO I=1,N
             CONT(I)=Y(I)+CONT(I)
          END DO
          CALL FCN(N,X,CONT,FF,RPAR,IPAR)
          NFCN=NFCN+1
          DO I=1,N
             CONT(I)=FF(I)+FF(I+N)
          END DO
          GOTO (31,32,31,32,31,32,33,55,55,55,41,42,41,42,41), IJOB
C ------ FULL MATRIX OPTION
 31      CONTINUE
         CALL SOL (N,LDE1,E1,CONT,IP1)
          GOTO 88
C ------ FULL MATRIX OPTION, SECOND ORDER
 41      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=1,NM1
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1)
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
          GOTO 88
C ------ BANDED MATRIX OPTION
 32      CONTINUE
         CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
          GOTO 88
C ------ BANDED MATRIX OPTION, SECOND ORDER
 42      CONTINUE
         DO J=1,M2
            SUM1=0.D0
            DO K=MM-1,0,-1
               SUM1=(CONT(J+K*M2)+SUM1)/FAC1
               DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
                  IM1=I+M1
                  CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
               END DO
            END DO
         END DO
         CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
         DO I=M1,1,-1
            CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
         END DO
          GOTO 88
C ------ HESSENBERG MATRIX OPTION
  33      CONTINUE
          DO MM=N-2,1,-1
             MP=N-MM
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 510
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 510         CONTINUE
             DO I=MP+1,N
                CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
             END DO
          END DO
          CALL SOLH(N,LDE1,E1,1,CONT,IP1)
          DO MM=1,N-2
             MP=N-MM
             DO I=MP+1,N
                CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
             END DO
             I=IPHES(MP)
             IF (I.EQ.MP) GOTO 640
             ZSAFE=CONT(MP)
             CONT(MP)=CONT(I)
             CONT(I)=ZSAFE
 640         CONTINUE
          END DO
C -----------------------------------
  88      CONTINUE
          ERR=0.D0
          DO I=1,N
             ERR=ERR+(CONT(I)/SCAL(I))**2
          END DO
          ERR=MAX(SQRT(ERR/N),1.D-10)
       END IF
       RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
       END
C
C     END OF SUBROUTINE ESTRAV
C
C ***********************************************************
C
      SUBROUTINE SLVROD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,E,LDE,IP,DY,AK,FX,YNEW,HD,IJOB,STAGE1)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E(LDE,NM1),
     &          IP(NM1),DY(N),AK(N),FX(N),YNEW(N)
      LOGICAL STAGE1
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      IF (HD.EQ.0.D0) THEN
         DO  I=1,N
           AK(I)=DY(I)
         END DO
      ELSE
         DO I=1,N
            AK(I)=DY(I)+HD*FX(I)
         END DO
      END IF
C
      GOTO (1,2,3,4,5,6,55,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
      CALL SOL (N,LDE,E,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
 48   MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(AK(JKM)+SUM)/FAC1
            DO I=1,NM1
               IM1=I+M1
               AK(IM1)=AK(IM1)+FJAC(I,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE,E,AK(M1+1),IP)
      DO I=M1,1,-1
         AK(I)=(AK(I)+AK(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
      CALL SOLB (N,LDE,E,MLE,MUE,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,N
            AK(I)=AK(I)+YNEW(I)
         END DO
      END IF
  45  MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(AK(JKM)+SUM)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               AK(IM1)=AK(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE,E,MLE,MUE,AK(M1+1),IP)
      DO I=M1,1,-1
         AK(I)=(AK(I)+AK(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      IF (STAGE1) THEN
      DO  I=1,N
         SUM=0.D0
         DO  J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
         END DO
         AK(I)=AK(I)+SUM
      END DO
      END IF
      CALL SOL (N,LDE,E,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,M1
            AK(I)=AK(I)+YNEW(I)
         END DO
         DO I=1,NM1
            SUM=0.D0
            DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
                SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J+M1)
            END DO
            IM1=I+M1
            AK(IM1)=AK(IM1)+SUM
         END DO
      END IF
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      IF (STAGE1) THEN
      DO I=1,N
         SUM=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*YNEW(J)
         END DO
         AK(I)=AK(I)+SUM
      END DO
      END IF
      CALL SOLB (N,LDE,E,MLE,MUE,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      IF (STAGE1) THEN
      DO I=1,N
         SUM=0.D0
         DO J=1,N
            SUM=SUM+FMAS(I,J)*YNEW(J)
         END DO
         AK(I)=AK(I)+SUM
      END DO
      END IF
      CALL SOL (N,LDE,E,AK,IP)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      IF (STAGE1) THEN
         DO I=1,M1
            AK(I)=AK(I)+YNEW(I)
         END DO
         DO I=1,NM1
            SUM=0.D0
            DO J=1,NM1
               SUM=SUM+FMAS(I,J)*YNEW(J+M1)
            END DO
            IM1=I+M1
            AK(IM1)=AK(IM1)+SUM
         END DO
      END IF
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      IF (STAGE1) THEN
      DO 624 I=1,N
         SUM=0.D0
         DO 623 J=1,N
  623       SUM=SUM+FMAS(I,J)*YNEW(J)
  624    AK(I)=AK(I)+SUM
      CALL SOLB (N,LDE,E,MLE,MUE,AK,IP)
      END IF
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVROD
C
C
C ***********************************************************
C
      SUBROUTINE SLVSEU(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,E,LDE,IP,IPHES,DEL,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FJAC(LDJAC,N),FMAS(LDMAS,NM1),E(LDE,NM1),DEL(N)
      DIMENSION IP(NM1),IPHES(N)
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,1,2,1,55,7,55,55,55,11,12,11,12,11), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      CALL SOL (N,LDE,E,DEL,IP)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(DEL(JKM)+SUM)/FAC1
            DO I=1,NM1
               IM1=I+M1
               DEL(IM1)=DEL(IM1)+FJAC(I,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE,E,DEL(M1+1),IP)
      DO I=M1,1,-1
         DEL(I)=(DEL(I)+DEL(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      CALL SOLB (N,LDE,E,MLE,MUE,DEL,IP)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      MM=M1/M2
      DO J=1,M2
         SUM=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM=(DEL(JKM)+SUM)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               DEL(IM1)=DEL(IM1)+FJAC(I+MUJAC+1-J,JKM)*SUM
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE,E,MLE,MUE,DEL(M1+1),IP)
      DO I=M1,1,-1
         DEL(I)=(DEL(I)+DEL(M2+I))/FAC1
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  HESSENBERG OPTION
      DO MMM=N-2,1,-1
         MP=N-MMM
         MP1=MP-1
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 110
         ZSAFE=DEL(MP)
         DEL(MP)=DEL(I)
         DEL(I)=ZSAFE
 110     CONTINUE
         DO I=MP+1,N
            DEL(I)=DEL(I)-FJAC(I,MP1)*DEL(MP)
         END DO
      END DO
      CALL SOLH(N,LDE,E,1,DEL,IP)
      DO MMM=1,N-2
         MP=N-MMM
         MP1=MP-1
         DO I=MP+1,N
            DEL(I)=DEL(I)+FJAC(I,MP1)*DEL(MP)
         END DO
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 240
         ZSAFE=DEL(MP)
         DEL(MP)=DEL(I)
         DEL(I)=ZSAFE
 240     CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVSEU
C
