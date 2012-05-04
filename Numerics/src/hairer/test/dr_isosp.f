        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NQDIM=200,NVDIM=200,NUDIM=1)
        PARAMETER(LWORK=500000,LIWORK=900000)
        DIMENSION Q(NQDIM),V(NVDIM),U(NUDIM),WORK(LWORK)
        DIMENSION IWORK(LIWORK),a(nvdim),rlam(68)
        EXTERNAL FPROB2,SOLO
           NUISO = 31
           NISO = NUISO + 1
           NBODY = NISO + 1
           NMRC = 1
           NBLK = NBODY*3 + 2
           NZGMAX = 8*NBODY + 2
           NZFMAX = 1
           CALL DATINS(NQ,NV,NU,NL,liwork,iwork,X0,Q,V,U)
C
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
           IWORK(23)=0       
           IWORK(24)=0 
           IWORK(1)=NZGMAX
           IWORK(2)=NZGMAX
           IWORK(14)=4
           IWORK(21)=NMRC
           IWORK(22)=NBLK
           IWORK(25)=100000
           IWORK(26)=100000
           IWORK(27)=2
           IWORK(28)=6
           WORK(5)=0.D0
           WORK(2)=0.9D0
           work(7)=0.d0
           work(8)=1/5.d0
           iwork(13)=0
         CALL HEM5(NQ,NV,NU,NL,FPROB2, 
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
c        WRITE(6,9921)(Q(I),I=1,NQ)
c        WRITE(6,9921)(V(I),I=1,NV) 
        WRITE(6,*)' STATISTICS : NFCN=',NFCN,' NSTEP=',NSTEP,' NACCPT=',
     &      NACCPT,' NREJCT=',NREJCT,' NDEC=',NDEC,' NSOL=',NSOL
 9921      FORMAT(1X,F22.16) 
        END
C
      SUBROUTINE SOLO(MODE,NR,NQ,NV,NU,NL,NZG,NZF,NZA,LRDO,
     &  LIDO,FPROB2,Q,V,U,DOWK,IDOWK)
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
     &              LRDO,LIDO,IDOWK,s,DOWK,FPROB2)
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
c
      SUBROUTINE FPROB2(IFCN,NQ,NV,NU,NL,NZG,NZF,lrda,NBLK,NMRC,
     &   NPGP,NPFL,INDGR,INDGC,INDFLR,INDFLC,
     &  T,P,V,U,XL,G,GP,F,GPP,GT,FL,QDOT,UDOT,AM)	
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER INDGR(NZG),INDGC(NZG),INDFLR(NZF),INDFLC(NZF)
      DIMENSION P(NQ),AM(lrda,Nmrc),V(NV),U(NU),UDOT(NU),QDOT(NQ),
     &     G(NL),GP(nzg),GPP(NL),F(NV),GT(NL),FL(nzf),XL(NL)
C
      COMMON / ISO / PMASS0,PMASS(100),PIN(100),PC(100),PD(100),
     $               PMASST,PINT,PCT,PDT,PAT,PE,PA,PRHO,
     $               PF0,PCL,PCQ,PFH,PXA,PYA,NUISO,NISO,NBODY
C  Dimensions
      NUISO = 31
c      NUISO = 31
      NISO = NUISO + 1
      NBODY = NISO + 1
      NP=NQ
C
      IF ((IFCN.EQ.1).OR.(IFCN.GE.7)) THEN
         DO 1010 J=1,NMRC
            DO 1000 I=1,NMRC*NBLK
               AM(I,J) = 0.D0
 1000       CONTINUE
 1010    CONTINUE
         AM(1,1) = PMASS0
         AM(2,1) = PMASS0
C
         DO 1020 I=1,NBODY
            L = 3*(I-1) + 2
            AM(L+1,1) = PMASS(I)
            AM(L+2,1) = PMASS(I)
            AM(L+3,1) = PIN(I)
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
C
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
       DO  J = 1,Nzg
        GP(J) = 0.D0
       end do
C
C  Rope
C
         GP(1) = -1.D0
         GP(2) = -1.D0
	 indgr(1)=1
         indgr(2)=2
         indgc(1)=1
         indgc(2)=2
         L=2
C
C  All bodies
C
         DO  I=1,NBODY
            L=L+1
            IROW = (I-1)*2+1
            ICOL = (I-1)*3+3
            GP(L ) =  1.D0
            INDGR(L)=IROW
            INDGC(L)=ICOL
            L=L+1
            GP(L ) =  1.D0
            INDGR(L)=IROW+1
            INDGC(L)=ICOL+1
            L=L+1
            GP(L ) =  -1.D0
            INDGR(L)=IROW+2
            INDGC(L)=ICOL
            L=L+1
            GP(L ) =  -1.D0
            INDGR(L)=IROW+3
            INDGC(L)=ICOL+1
         end do
C
C  Triangle
C
         L=L+1
         GP(L) = PD(1)*COS(P(5))
         INDGR(L)=1
         INDGC(L)=5
         L=L+1
         GP(L) = PD(1)*SIN(P(5))
         INDGR(L)=2
         INDGC(L)=5
         L=L+1
         GP(L) = PD(1)*
     $        (SQRT(3.D0)*0.5D0*SIN(P(5)) + 0.5D0*COS(P(5)))
         INDGR(L)=3
         INDGC(L)=5
         L=L+1
         GP(L) = -PD(1)*
     $        (SQRT(3.D0)*0.5D0*COS(P(5)) - 0.5D0*SIN(P(5)))
         INDGR(L)=4
         INDGC(L)=5
C
C  Isolators
C
         DO  I=2,NBODY
            IROW = (I-1)*2 + 1
            ICOL = (I-1)*3 + 5
            L=L+1
            GP(L) = PC(I)*COS(P(ICOL))
            INDGR(L)=IROW
            INDGC(L)=ICOL
            L=L+1
            GP(L) = PC(I)*SIN(P(ICOL))
            INDGR(L)=IROW+1
            INDGC(L)=ICOL
            L=L+1
            GP(L) = PD(I)*COS(P(ICOL))
            INDGR(L)=IROW+2
            INDGC(L)=ICOL
            L=L+1
            GP(L) = PD(I)*SIN(P(ICOL))
            INDGR(L)=IROW+3
            INDGC(L)=ICOL
         END DO
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
C
