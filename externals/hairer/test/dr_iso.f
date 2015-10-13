C         include 'linsp.f'
C         include 'hem5.f'
C         include 'decsol.f'
C         include 'lapack.f'
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NQDIM=200,NVDIM=200,NUDIM=1)
        PARAMETER(LWORK=200000,LIWORK=9000)
        DIMENSION Q(NQDIM),V(NVDIM),U(NUDIM),WORK(LWORK)
        DIMENSION IWORK(LIWORK),a(nvdim),rlam(1)        
	EXTERNAL FPROB1,SOLO
        write(6,*)'LAPACK VERSION STANDARD'
         NUISO = 31
         NISO = NUISO + 1
         NBODY = NISO + 1
         NP = 2 + NBODY*3
         NV = NP
         NG = 2 + NBODY*2
         NL = NG
         NU=0
           CALL DATINS(NQ,NV,NU,NL,liwork,iwork,X0,Q,V,U)
           RTOL=1.d-5
           ATOL = 1.D-1*RTOL  
           ITOL=0   
           IOUT=0
           h0=1.d-3
           X = 0.D0
           XEND = 0.1d0
           H = H0  
           DO 10 I=1,30
           WORK(I)=0.D0
  10       IWORK(I)=0
           WORK(5)=0.D0
           WORK(2)=0.9D0
           work(7)=0.d0
           work(8)=1/5.d0
           iwork(14)=2
         CALL HEM5(NQ,NV,NU,NL,FPROB1, 
     &                X,Q,V,U,a,rlam,Xend,H,
     &                RTOL,ATOL,ITOL,
     &                SOLO,IOUT,
     &                WORK,LWORK,IWORK,LIWORK,IDID)
        NSTEP=IWORK(31)
	NACCPT=IWORK(32)
	NREJCT=IWORK(33)
	NFCN=IWORK(34) 
	NDEC=IWORK(37) 
        NSOL=IWORK(38)
        WRITE(6,9921)(Q(I),I=1,NQ)
        WRITE(6,9921)(V(I),I=1,NV) 
        WRITE(6,*)' STATISTICS : NFCN=',NFCN,' NSTEP=',NSTEP,' NACCPT=',
     &      NACCPT,' NREJCT=',NREJCT,' NDEC=',NDEC,' NSOL=',NSOL
 9921      FORMAT(1X,F22.16) 
        END
C
      SUBROUTINE SOLO(MODE,NR,NQ,NV,NU,NL,NZG,NZF,NZA,LRDO,
     &  LIDO,FPROB,Q,V,U,DOWK,IDOWK)
C --- BY USING "CONTD8", THE CONTINUOUS COLLOCATION SOLUTION
        IMPLICIT REAL*8 (A-H,O-Z) 
        parameter (nrpts=4)
        PARAMETER(LWORK=10000,LIWORK=1000)
        DIMENSION Q(NQ),V(NV),U(NU)
        DIMENSION DOWK(LRDO),IDOWK(LIDO)
        EXTERNAL Fprob
        XOLD=DOWK(1)
        H=DOWK(2) 
        ides=8             
           hstep=H/(nrpts-1)
          if (nr.gt.1) then   
           do  i=2,nrpts
           s=xold+(i-1)*hstep 
           tt=(s-xold)/h  
             xend=s
        write(6,*)'nr,s',nr,s
         v1= ct4(MODE,IDES,NQ,NV,NU,NL,NZG,NZF,NZA,
     &              LRDO,LIDO,IDOWK,s,DOWK,FPROB)
        write(6,*) 'v1',v1
	end do
         end if
        RETURN
        END
c
      SUBROUTINE DATINS(NP,NV,NU,NL,LILA,ILA,T,P,V,U) 
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ILA(LILA)
      DIMENSION P(NP),V(NV),U(NU)
C
      COMMON / ISO / PMASS0,PMASS(100),PIN(100),PC(100),PD(100),
     $               PMASST,PINT,PCT,PDT,PAT,PE,PA,PRHO,
     $               PF0,PCL,PCQ,PFH,PXA,PYA,NUISO,NISO,NBODY
C
C
C  Dimensions
C
C
      NUISO = 31
c      NUISO = 31
      NISO = NUISO + 1
      NBODY = NISO + 1
      NP = 2 + NBODY*3
      NV = NP
      NG = 2 + NBODY*2
      NL = NG
      NU=0
C Physical parameters
C
      PMASS0 = 0.D0
      PMASST = 34.D0
      PMASSI = 15.D0
      PMASSH = 9.8D0
      PINT = 3.1D0
      PINI = 0.35D0
      PINH = 0.05D0
      PCT = 0.37D0
      PDT = PCT
      PCI = 0.08D0
      PDI = 0.12D0
      PCH = 0.08D0
      PDH = 0.16D0
C
      PE = 8.D9
      PA = 3.4D-4
      PRHO = 3325.D0
      PF0 = 2.3D5
      PCL = SQRT(PE/PRHO)
      PCQ = SQRT(PF0/(PA*PRHO))
      PFH = PE*PA/PCL
      PAT = SQRT(3.D0)*PCT
      PXA = 0.D0
      PYA = 0.D0
C
      PMASS(1) = PMASST
      PIN(1) = PINT
      PC(1) = PCT
      PD(1) = PDT
C
      DO 1020 I=2,NISO
         PMASS(I) = PMASSI
         PIN(I) = PINI
         PC(I) = PCI
         PD(I) = PDI
 1020 CONTINUE
C
      PMASS(NBODY) = PMASSH
      PIN(NBODY) = PINH
      PC(NBODY) = PCH
      PD(NBODY) = PDH
C
C
C  Initial values
C
C   Position of isolators
C
      DO IPOS=1,NP
        P(IPOS)=0.D0
      END DO
      P(NP) = 0.0D0
      P(NP-1) = PYA - PD(NBODY)
      P(NP-2) = PXA
      DO 1030 L=NISO-1,1,-1
         IX = (L-1)*3 + 6
         IY = IX + 1
         IPHI = IX + 2
         P(IX) = PXA
         P(IY) = P(IY+3) - PC(L+2) - PD(L+1)
 1030 CONTINUE
C
C   Position of triangle
C
      P(3) = PXA - PAT/2.D0
      P(4) = P(7) - PC(2) - PCT/2.D0
      P(5) = 0.0D0
C
C   Position of rope
C
      P(1) = P(3)
      P(2) = P(4) - PCT
C
C   Velocities
      DO 1040 I=1,NV
         V(I) = 0.D0
 1040 CONTINUE
C
      RETURN
      END
c**********************************************************************
      SUBROUTINE FPROB1(IFCN,NQ,NV,NU,NL,NZG,NZF,lrda,NBLK,NMRC,
     &   NPGP,NPFL,INDGR,INDGC,INDFLR,INDFLC,
     &  T,P,V,U,XL,G,GP,F,GPP,GT,FL,QDOT,UDOT,AM)	
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER INDGR(NZG),INDGC(NZG),INDFLR(NZF),INDFLC(NZF)
      DIMENSION P(NQ),AM(lrda,NV),V(NV),U(NU),UDOT(NU),QDOT(NQ),
     &     G(NL),GP(NL,NV),GPP(NL),F(NV),GT(NL),FL(NV,NL),XL(NL)
C
      COMMON / ISO / PMASS0,PMASS(100),PIN(100),PC(100),PD(100),
     $               PMASST,PINT,PCT,PDT,PAT,PE,PA,PRHO,
     $               PF0,PCL,PCQ,PFH,PXA,PYA,NUISO,NISO,NBODY
C
C  Dimensions
      NUISO = 31
      NISO = NUISO + 1
      NBODY = NISO + 1
C
      IF ((IFCN.EQ.1).OR.(IFCN.GE.7)) THEN
        DO 10 I = 1,NV
        DO 10 J = 1,NV
  10      AM(I,J) = 0.D0 
         AM(1,1) = PMASS0
         AM(2,2) = PMASS0
C
         DO 1020 I=1,NBODY
            L = 3*(I-1) + 2
            AM(L+1,L+1) = PMASS(I)
            AM(L+2,L+2) = PMASS(I)
            AM(L+3,L+3) = PIN(I)
 1020    CONTINUE
      END IF
C
      IF ((IFCN.EQ.1).OR.(IFCN.EQ.5).OR.(IFCN.EQ.7).OR.(IFCN.EQ.8)) THEN
C
         BETA = -ATAN(V(1)/PCQ)
         FAC = PF0 + V(2)*PFH
         F(1) = FAC*SIN(BETA)
         F(2) = -FAC*COS(BETA)
         DO 1070 I=3,NV
            F(I) = 0.D0
 1070    CONTINUE
      END IF
C
      IF (IFCN.EQ.4) THEN
C  Seil zum Abstandshalter
C
         G(1) = P(3) + PC(1)*SIN(P(5)) - P(1)
         G(2) = P(4) - PC(1)*COS(P(5)) - P(2)
C
C  Abstandshalter zum Isolator
C
         G(3) = P(6) + PC(2)*SIN(P(8)) - P(3) - PD(1)*
     $        (SQRT(3.D0)*0.5D0*COS(P(5)) - 0.5D0*SIN(P(5)))
         G(4) = P(7) - PC(2)*COS(P(8)) - P(4) - PD(1)*
     $        (SQRT(3.D0)*0.5D0*SIN(P(5)) + 0.5D0*COS(P(5)))
C
C  Isolator-Gelenke
C
         DO 1090 I=2,NBODY-1
            JG = 2*(I-1) + 3
            LI = 3*(I-1) + 3
            LI1 = LI + 3
            G(JG)   = P(LI1) + PC(I+1)*SIN(P(LI1+2))
     $           -P(LI) + PD(I)*SIN(P(LI+2))
            G(JG+1) = P(LI1+1) - PC(I+1)*COS(P(LI1+2))
     $           -P(LI+1) - PD(I)*COS(P(LI+2))
 1090    CONTINUE
C
C  Aufhaengepunkt
C
         JG = 2*(NBODY-1) + 3
         LI = 3*(NBODY-1) + 3
         LI1 = LI + 3
         G(JG)   = -P(LI) + PD(NBODY)*SIN(P(LI+2))
         G(JG+1) = -P(LI+1) - PD(NBODY)*COS(P(LI+2))
      END IF
C      
      IF ((IFCN.EQ.6).OR.(IFCN.GE.10)) THEN
       DO  J = 1,NV
        DO I=1,NL
        GP(I,J) = 0.D0
       end do
       END DO
C
C  Rope
C
         GP(1,1) = -1.D0
         GP(2,2) = -1.D0
C
C  All bodies
C
         DO 1050 I=1,NBODY
            IROW = (I-1)*2+1
            ICOL = (I-1)*3+3
            GP(IROW  , ICOL  ) =  1.D0
            GP(IROW+1, ICOL+1) =  1.D0
            GP(IROW+2, ICOL  ) = -1.D0
            GP(IROW+3, ICOL+1) = -1.D0
 1050    CONTINUE
C
C  Triangle
C
         GP(1,5) = PD(1)*COS(P(5))
         GP(2,5) = PD(1)*SIN(P(5))
         GP(3,5) = PD(1)*
     $        (SQRT(3.D0)*0.5D0*SIN(P(5)) + 0.5D0*COS(P(5)))
         GP(4,5) = -PD(1)*
     $        (SQRT(3.D0)*0.5D0*COS(P(5)) - 0.5D0*SIN(P(5)))
C
C  Isolators
C
         DO 1060 I=2,NBODY
            IROW = (I-1)*2 + 1
            ICOL = (I-1)*3 + 5
C           WRITE(6,*) ' IROW, ICOL', IROW, ICOL
C           WRITE(6,*) ' NP, NG', NP, NG
            GP(IROW  , ICOL) = PC(I)*COS(P(ICOL))
            GP(IROW+1, ICOL) = PC(I)*SIN(P(ICOL))
            GP(IROW+2, ICOL) = PD(I)*COS(P(ICOL))
            GP(IROW+3, ICOL) = PD(I)*SIN(P(ICOL))
 1060    CONTINUE
      END IF
C
      IF ((IFCN.EQ.5).OR.(IFCN.EQ.7)) THEN
	WRITE(6,*)'GPP NOT AVAILABLE'
      END IF
C
      IF ((IFCN.EQ.3).OR.(IFCN.EQ.6).OR.(IFCN.GE.10)) THEN
      DO  I=1,NL
         GT(I)=0.
      END DO
      END IF
C
      IF (IFCN.EQ.0) THEN
      DO  I=1,NU
         UDOT(I)=0.
      END DO
      END IF
C
      IF ((IFCN.EQ.1).OR.(IFCN.EQ.2).OR.(IFCN.EQ.10)) THEN
      DO  I=1,NQ
         QDOT(I)=V(I)
      END DO
      END IF
      RETURN
      END

