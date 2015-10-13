C The codes MA28**, MA30**, MC*** are part of the Harwell Subroutine Library,
C which is developed
C and maintained by staff at the Harwell Laboratory of the United Kingdom
C Atomic Energy Authority. The researcher may use this code SOLELY for
C research purposes subject to the conditions below.
C 
C 
C  On receipt of this code the requester must return the following 
C form (signed) to Mr. Sid Marlow, Computer Science and Systems Division,
C B7.12, AERE Harwell, Didcot, Oxon, England OX11 0RA.
C 
C 
C The conditions as specified below are accepted and shall be complied with.
C 
C 
C Signature ...........................   Date ............................
C 
C Name (BLOCK CAPITALS please) ............................................
C 
C on behalf of (BLOCK CAPITALS please) ....................................
C 
C Address..................................................................
C 
C  ........................................................................
C 
C Position ................................................................
C 
C  
C  The conditions of use of the Harwell Subroutine Library or any part
C thereof by any external user are as follows:
C 
C (i)   The subroutines may only be used for research purposes by the person or
C       organization to whom they are supplied ("the User"), and must not be
C       copied by the User for use by any other persons or organizations.
C       All information on the subroutines is  provided by the Atomic Energy
C       Authority ("the Authority") to the User on the understanding that
C       the details thereof are confidential.
C 
C (ii)  All publications issued by the User which include results obtained
C       with the help of one or more of the subroutines supplied shall
C       acknowledge the use of the subroutines.
C 
C (iii) The subroutines supplied may be modified by or on behalf of the User
C       for such use in research applications but at no time shall such
C       subroutines or modifications thereof become the property of the User.
C       The Authority shall not be liable for any direct or consequential
C       loss or damage whatsoever arising out of the use of subroutines by
C       the User.
C 
C (iv)  Any use of the subroutines supplied in any commercial application
C       shall be subject to prior written agreement between the Authority
C       and the User on suitable terms and conditions (including financial
C       conditions).
C 
c ---------------------------------------------------------
C
C########################################################################
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias ma28ad ma28bd ma28cd
c###### calls   ma30    mc20    mc22    mc23    mc24
      SUBROUTINE MA28AD(N, NZ, A, LICN, IRN, LIRN, ICN, U, IKEEP, IW, W,
     * IFLAG)
c this subroutine performs the lu factorization of a.
c
c the parameters are as follows.....
c n     order of matrix  not altered by subroutine.
c nz    number of non-zeros in input matrix  not altered by subroutine.
c a is a  real array  length licn.  holds non-zeros of matrix on entry
c     and non-zeros of factors on exit.  reordered by mc20a/ad and
c     mc23a/ad and altered by ma30a/ad.
c licn  integer  length of arrays a and icn.  not altered by subroutine.
c irn   integer array of length lirn.  holds row indices on input.
c     used as workspace by ma30a/ad to hold column orientation of
c     matrix.
c lirn  integer  length of array irn. not altered by the subroutine.
c icn   integer array of length licn.  holds column indices on entry
c     and column indices of decomposed matrix on exit. reordered by
c     mc20a/ad and mc23a/ad and altered by ma30a/ad.
c u     real variable  set by user to control bias towards numeric or
c     sparsity pivoting.  u=1.0 gives partial pivoting while u=0. does
c     not check multipliers at all.  values of u greater than one are
c     treated as one while negative values are treated as zero.  not
c     altered by subroutine.
c ikeep  integer array of length 5*n  used as workspace by ma28a/ad
c     (see later comments).  it is not required to be set on entry
c     and, on exit, it contains information about the decomposition.
c     it should be preserved between this call and subsequent calls
c     to ma28b/bd or ma28c/cd.
c     ikeep(i,1),i=1,n  holds the total length of the part of row i
c     in the diagonal block.
c     row ikeep(i,2),i=1,n  of the input matrix is the ith row in
c     pivot order.
c     column ikeep(i,3),i=1,n  of the input matrix is the ith column
c     in pivot order.
c     ikeep(i,4),i=1,n  holds the length of the part of row i in
c     the l part of the l/u decomposition.
c     ikeep(i,5),i=1,n  holds the length of the part of row i in the
c     off-diagonal blocks.  if there is only one diagonal block,
c     ikeep(1,5) will be set to -1.
c iw    integer array of length 8*n.  if the option nsrch.le.n is
c     used, then the length of array iw can be reduced to 7*n.
c w real array  length n.  used by mc24a/ad both as workspace and to
c     return growth estimate in w(1).  the use of this array by ma28a/ad
c     is thus optional depending on common block logical variable grow.
c iflag  integer variable  used as error flag by routine.  a positive
c     or zero value on exit indicates success.  possible negative
c     values are -1 through -14.
c
      INTEGER N, NZ, LICN, LIRN, IFLAG
      INTEGER IRN(LIRN), ICN(LICN), IKEEP(N,5), IW(N,8)
      DOUBLE PRECISION A(LICN), U, W(N)
c
c common and private variables.
c     common block ma28f/fd is used merely
c     to communicate with common block ma30f/fd  so that the user
c     need not declare this common block in his main program.
c the common block variables are as follows ...
c lp,mp  integer  default value 6 (line printer).  unit number
c     for error messages and duplicate element warning resp.
c nlp,mlp  integer  unit number for messages from ma30a/ad and
c     mc23a/ad resp.  set by ma28a/ad to value of lp.
c lblock  logical  default value true.  if true mc23a/ad is used
c     to first permute the matrix to block lower triangular form.
c grow    logical  default value true.  if true then an estimate
c     of the increase in size of matrix elements during l/u
c     decomposition is given by mc24a/ad.
c eps,rmin,resid  real/double precision variables not referenced
c     by ma28a/ad.
c irncp,icncp  integer  set to number of compresses on arrays irn and
c     icn/a respectively.
c minirn,minicn  integer  minimum length of arrays irn and icn/a
c     respectively, for success on future runs.
c irank  integer   estimated rank of matrix.
c mirncp,micncp,mirank,mirn,micn integer variables.  used to
c     communicate between ma30f/fd and ma28f/fd values of abovenamed
c     variables with somewhat similar names.
c abort1,abort2  logical variables with default value true.  if false
c     then decomposition will be performed even if the matrix is
c     structurally or numerically singular respectively.
c aborta,abortb  logical variables used to communicate values of
c     abort1 and abort2 to ma30a/ad.
c abort  logical  used to communicate value of abort1 to mc23a/ad.
c abort3  logical variable not referenced by ma28a/ad.
c idisp   integer array  length 2.  used to communicate information
c     on decomposition between this call to ma28a/ad and subsequent
c     calls to ma28b/bd and ma28c/cd.  on exit, idisp(1) and
c     idisp(2) indicate position in arrays a and icn of the
c     first and last elements in the l/u decomposition of the
c     diagonal blocks, respectively.
c numnz  integer  structural rank of matrix.
c num    integer  number of diagonal blocks.
c large  integer  size of largest diagonal block.
c
c see block data for further comments on common block variables.
c see code for comments on private variables.
c
      DOUBLE PRECISION TOL, THEMAX, BIG, DXMAX, ERRMAX, DRES, CGCE,
     * TOL1, BIG1, UPRIV, RMIN, EPS, RESID, ZERO
      INTEGER IDISP(2)
      LOGICAL GROW, LBLOCK, ABORT, ABORT1, ABORT2, ABORT3, ABORTA,
     * ABORTB, LBIG, LBIG1
      COMMON /MA28ED/ LP, MP, LBLOCK, GROW
      COMMON /MA28FD/ EPS, RMIN, RESID, IRNCP, ICNCP, MINIRN, MINICN,
     * IRANK, ABORT1, ABORT2
      COMMON /MA28GD/ IDISP
      COMMON /MA28HD/ TOL, THEMAX, BIG, DXMAX, ERRMAX, DRES, CGCE,
     * NDROP, MAXIT, NOITER, NSRCH, ISTART, LBIG
      COMMON /MA30ID/ TOL1, BIG1, NDROP1, NSRCH1, LBIG1
      COMMON /MA30ED/ NLP, ABORTA, ABORTB, ABORT3
      COMMON /MA30FD/ MIRNCP, MICNCP, MIRANK, MIRN, MICN
      COMMON /MC23BD/ MLP, NUMNZ, NUM, LARGE, ABORT
      COMMON /LPIVOT/ LPIV(10),LNPIV(10),MAPIV,MANPIV,IAVPIV,
     *                IANPIV,KOUNTL
c
c some  initialization and transfer of information between
c     common blocks (see earlier comments).
      DATA ZERO /0.0D0/
      IFLAG = 0
      ABORTA = ABORT1
      ABORTB = ABORT2
      ABORT = ABORT1
      MLP = LP
      NLP = LP
      TOL1 = TOL
      LBIG1 = LBIG
      NSRCH1 = NSRCH
c upriv private copy of u is used in case it is outside
c     range  zero to one  and  is thus altered by ma30a/ad.
      UPRIV = U
c simple data check on input variables and array dimensions.
      IF (N.GT.0) GO TO 10
      IFLAG = -8
      IF (LP.NE.0) WRITE (LP,99999) N
      GO TO 210
   10 IF (NZ.GT.0) GO TO 20
      IFLAG = -9
      IF (LP.NE.0) WRITE (LP,99998) NZ
      GO TO 210
   20 IF (LICN.GE.NZ) GO TO 30
      IFLAG = -10
      IF (LP.NE.0) WRITE (LP,99997) LICN
      GO TO 210
   30 IF (LIRN.GE.NZ) GO TO 40
      IFLAG = -11
      IF (LP.NE.0) WRITE (LP,99996) LIRN
      GO TO 210
c
c data check to see if all indices lie between 1 and n.
   40 DO 50 I=1,NZ
        IF (IRN(I).GT.0 .AND. IRN(I).LE.N .AND. ICN(I).GT.0 .AND.
     *   ICN(I).LE.N) GO TO 50
        IF (IFLAG.EQ.0 .AND. LP.NE.0) WRITE (LP,99995)
        IFLAG = -12
        IF (LP.NE.0) WRITE (LP,99994) I, A(I), IRN(I), ICN(I)
   50 CONTINUE
      IF (IFLAG.LT.0) GO TO 220
c
c sort matrix into row order.
      CALL MC20AD(N, NZ, A, ICN, IW, IRN, 0)
c part of ikeep is used here as a work-array.  ikeep(i,2) is
c     the last row to have a non-zero in column i.  ikeep(i,3)
c     is the off-set of column i from the start of the row.
      DO 60 I=1,N
        IKEEP(I,2) = 0
        IKEEP(I,1) = 0
   60 CONTINUE
c
c check for duplicate elements .. summing any such entries and
c     printing a warning message on unit mp.
c move is equal to the number of duplicate elements found.
      MOVE = 0
c the loop also calculates the largest element in the matrix, themax.
      THEMAX = ZERO
c j1 is position in arrays of first non-zero in row.
      J1 = IW(1,1)
      DO 130 I=1,N
        IEND = NZ + 1
        IF (I.NE.N) IEND = IW(I+1,1)
        LENGTH = IEND - J1
        IF (LENGTH.EQ.0) GO TO 130
        J2 = IEND - 1
        NEWJ1 = J1 - MOVE
        DO 120 JJ=J1,J2
          J = ICN(JJ)
          THEMAX = DMAX1(THEMAX,DABS(A(JJ)))
          IF (IKEEP(J,2).EQ.I) GO TO 110
c first time column has ocurred in current row.
          IKEEP(J,2) = I
          IKEEP(J,3) = JJ - MOVE - NEWJ1
          IF (MOVE.EQ.0) GO TO 120
c shift necessary because of  previous duplicate element.
          NEWPOS = JJ - MOVE
          A(NEWPOS) = A(JJ)
          ICN(NEWPOS) = ICN(JJ)
          GO TO 120
c duplicate element.
  110     MOVE = MOVE + 1
          LENGTH = LENGTH - 1
          JAY = IKEEP(J,3) + NEWJ1
          IF (MP.NE.0) WRITE (MP,99993) I, J, A(JJ)
          A(JAY) = A(JAY) + A(JJ)
          THEMAX = DMAX1(THEMAX,DABS(A(JAY)))
  120   CONTINUE
        IKEEP(I,1) = LENGTH
        J1 = IEND
  130 CONTINUE
c
c knum is actual number of non-zeros in matrix with any multiple
c     entries counted only once.
      KNUM = NZ - MOVE
      IF (.NOT.LBLOCK) GO TO 140
c
c perform block triangularisation.
      CALL MC23AD(N, ICN, A, LICN, IKEEP, IDISP, IKEEP(1,2),
     *IKEEP(1,3), IKEEP(1,5), IW(1,3), IW)
      IF (IDISP(1).GT.0) GO TO 170
      IFLAG = -7
      IF (IDISP(1).EQ.-1) IFLAG = -1
      IF (LP.NE.0) WRITE (LP,99992)
      GO TO 210
c
c block triangularization not requested.
c move structure to end of data arrays in preparation for
c     ma30a/ad.
c also set lenoff(1) to -1 and set permutation arrays.
  140 DO 150 I=1,KNUM
        II = KNUM - I + 1
        NEWPOS = LICN - I + 1
        ICN(NEWPOS) = ICN(II)
        A(NEWPOS) = A(II)
  150 CONTINUE
      IDISP(1) = 1
      IDISP(2) = LICN - KNUM + 1
      DO 160 I=1,N
        IKEEP(I,2) = I
        IKEEP(I,3) = I
  160 CONTINUE
      IKEEP(1,5) = -1
  170 IF (LBIG) BIG1 = THEMAX
      IF (NSRCH.LE.N) GO TO 180
c
c perform l/u decomosition on diagonal blocks.
      CALL MA30AD(N, ICN, A, LICN, IKEEP, IKEEP(1,4), IDISP,
     *IKEEP(1,2), IKEEP(1,3), IRN, LIRN, IW(1,2), IW(1,3), IW(1,4),
     *IW(1,5), IW(1,6), IW(1,7), IW(1,8), IW, UPRIV, IFLAG)
      GO TO 190
c this call if used if nsrch has been set less than or equal n.
c     in this case, two integer work arrays of length can be saved.
  180 CALL MA30AD(N, ICN, A, LICN, IKEEP, IKEEP(1,4), IDISP,
     * IKEEP(1,2), IKEEP(1,3), IRN, LIRN, IW(1,2), IW(1,3), IW(1,4),
     * IW(1,5), IW, IW, IW(1,6), IW, UPRIV, IFLAG)
c
c transfer common block information.
  190 MINIRN = MAX0(MIRN,NZ)
      MINICN = MAX0(MICN,NZ)
      IRNCP = MIRNCP
      ICNCP = MICNCP
      IRANK = MIRANK
      NDROP = NDROP1
      IF (LBIG) BIG = BIG1
      IF (IFLAG.GE.0) GO TO 200
      IF (LP.NE.0) WRITE (LP,99991)
      GO TO 210
c
c reorder off-diagonal blocks according to pivot permutation.
  200 I1 = IDISP(1) - 1
      IF (I1.NE.0) CALL MC22AD(N, ICN, A, I1, IKEEP(1,5), IKEEP(1,2),
     * IKEEP(1,3), IW, IRN)
      I1 = IDISP(1)
      IEND = LICN - I1 + 1
c
c optionally calculate element growth estimate.
      IF (GROW) CALL MC24AD(N, ICN, A(I1), IEND, IKEEP, IKEEP(1,4), W)
c increment growth estimate by original maximum element.
      IF (GROW) W(1) = W(1) + THEMAX
      IF (GROW .AND. N.GT.1) W(2) = THEMAX
c set flag if the only error is due to duplicate elements.
      IF (IFLAG.GE.0 .AND. MOVE.NE.0) IFLAG = -14
      GO TO 220
  210 IF (LP.NE.0) WRITE (LP,99990)
  220 RETURN
99999 FORMAT (36X, 17HN OUT OF RANGE = , I10)
99998 FORMAT (36X, 18HNZ NON POSITIVE = , I10)
99997 FORMAT (36X, 17HLICN TOO SMALL = , I10)
99996 FORMAT (36X, 17HLIRN TOO SMALL = , I10)
99995 FORMAT (54H ERROR RETURN FROM MA28A/AD BECAUSE INDICES FOUND OUT ,
     * 8HOF RANGE)
99994 FORMAT (1X, I6, 22HTH ELEMENT WITH VALUE , 1PD22.14, 9H IS OUT O,
     * 21HF RANGE WITH INDICES , I8, 2H ,, I8)
99993 FORMAT (31H DUPLICATE ELEMENT IN POSITION , I8, 2H ,, I8,
     * 12H WITH VALUE , 1PD22.14)
99992 FORMAT (36X, 26HERROR RETURN FROM MC23A/AD)
99991 FORMAT (36X, 26HERROR RETURN FROM MA30A/AD)
99990 FORMAT (36H+ERROR RETURN FROM MA28A/AD BECAUSE )
      END
c
c==============================================================
c
      SUBROUTINE MA28BD(N, NZ, A, LICN, IVECT, JVECT, ICN, IKEEP, IW, W,
     * IFLAG)
c this subroutine factorizes a matrix of a similar sparsity
c     pattern to that previously factorized by ma28a/ad.
c the parameters are as follows ...
c n      integer  order of matrix  not altered by subroutine.
c nz     integer  number of non-zeros in input matrix  not altered
c     by subroutine.
c a      real/double precision array  length licn.  holds non-zeros of
c     matrix on entry and non-zeros of factors on exit.  reordered by
c     ma28d/dd and altered by subroutine ma30b/bd.
c licn   integer  length of arrays a and icn.  not altered by
c     subroutine.
c ivect,jvect  integer arrays of length nz.  hold row and column
c     indices of non-zeros respectively.  not altered by subroutine.
c icn    integer array of length licn.  same array as output from
c     ma28a/ad.  unchanged by ma28b/bd.
c ikeep  integer array of length 5*n.  same array as output from
c     ma28a/ad.  unchanged by ma28b/bd.
c iw     integer array  length 5*n.  used as workspace by ma28d/dd and
c     ma30b/bd.
c w      real/double precision array  length n.  used as workspace
c     by ma28d/dd,ma30b/bd and (optionally) mc24a/ad.
c iflag  integer  used as error flag with positive or zero value
c     indicating success.
c
      INTEGER N, NZ, LICN, IW(N,5), IFLAG
      INTEGER IKEEP(N,5), IVECT(NZ), JVECT(NZ), ICN(LICN)
      DOUBLE PRECISION A(LICN), W(N)
c
c private and common variables.
c unless otherwise stated common block variables are as in ma28a/ad.
c     those variables referenced by ma28b/bd are mentioned below.
c lp,mp  integers  used as in ma28a/ad as unit number for error and
c     warning messages, respectively.
c nlp    integer variable used to give value of lp to ma30e/ed.
c eps    real/double precision  ma30b/bd will output a positive value
c     for iflag if any modulus of the ratio of pivot element to the
c     largest element in its row (u part only) is less than eps (unless
c     eps is greater than 1.0 when no action takes place).
c rmin   real/double precision  variable equal to the value of this
c     minimum ratio in cases where eps is less than or equal to 1.0.
c meps,mrmin  real/double precision variables used by the subroutine
c     to communicate between common blocks ma28f/fd and ma30g/gd.
c idisp  integer array  length 2  the same as that used by ma28a/ad.
c     it is unchanged by ma28b/bd.
c
c see block data or ma28a/ad for further comments on variables
c     in common.
c see code for comments on private variables.
c
      LOGICAL GROW, LBLOCK, ABORTA, ABORTB, ABORT1, ABORT2, ABORT3,
     * LBIG, LBIG1
      INTEGER IDISP(2)
      DOUBLE PRECISION EPS, MEPS, RMIN, MRMIN, RESID, TOL,
     * THEMAX, BIG, DXMAX, ERRMAX, DRES, CGCE, TOL1, BIG1
c
      COMMON /MA28ED/ MP, LP, LBLOCK, GROW
      COMMON /MA28FD/ EPS, RMIN, RESID, IRNCP, ICNCP, MINIRN, MINICN,
     * IRANK, ABORT1, ABORT2
      COMMON /MA28GD/ IDISP
      COMMON /MA28HD/ TOL, THEMAX, BIG, DXMAX, ERRMAX, DRES, CGCE,
     * NDROP, MAXIT, NOITER, NSRCH, ISTART, LBIG
      COMMON /MA30ED/ NLP, ABORTA, ABORTB, ABORT3
      COMMON /MA30GD/ MEPS, MRMIN
      COMMON /MA30ID/ TOL1, BIG1, NDROP1, NSRCH1, LBIG1
c
c check to see if elements were dropped in previous ma28a/ad call.
      IF (NDROP.EQ.0) GO TO 10
      IFLAG = -15
      WRITE (6,99999) IFLAG, NDROP
      GO TO 70
   10 IFLAG = 0
      MEPS = EPS
      NLP = LP
c simple data check on variables.
      IF (N.GT.0) GO TO 20
      IFLAG = -11
      IF (LP.NE.0) WRITE (LP,99998) N
      GO TO 60
   20 IF (NZ.GT.0) GO TO 30
      IFLAG = -10
      IF (LP.NE.0) WRITE (LP,99997) NZ
      GO TO 60
   30 IF (LICN.GE.NZ) GO TO 40
      IFLAG = -9
      IF (LP.NE.0) WRITE (LP,99996) LICN
      GO TO 60
c
   40 CALL MA28DD(N, A, LICN, IVECT, JVECT, NZ, ICN, IKEEP, IKEEP(1,4),
     * IKEEP(1,5), IKEEP(1,2), IKEEP(1,3), IW(1,3), IW, W(1), IFLAG)
c themax is largest element in matrix.
      THEMAX = W(1)
      IF (LBIG) BIG1 = THEMAX
c idup equals one if there were duplicate elements, zero otherwise.
      IDUP = 0
      IF (IFLAG.EQ.(N+1)) IDUP = 1
      IF (IFLAG.LT.0) GO TO 60
c
c perform row-gauss elimination on the structure received from ma28d/dd
      CALL MA30BD(N, ICN, A, LICN, IKEEP, IKEEP(1,4), IDISP,
     * IKEEP(1,2), IKEEP(1,3), W, IW, IFLAG)
c
c transfer common block information.
      IF (LBIG) BIG1 = BIG
      RMIN = MRMIN
      IF (IFLAG.GE.0) GO TO 50
      IFLAG = -2
      IF (LP.NE.0) WRITE (LP,99995)
      GO TO 60
c
c optionally calculate the growth parameter.
   50 I1 = IDISP(1)
      IEND = LICN - I1 + 1
      IF (GROW) CALL MC24AD(N, ICN, A(I1), IEND, IKEEP, IKEEP(1,4), W)
c increment estimate by largest element in input matrix.
      IF (GROW) W(1) = W(1) + THEMAX
      IF (GROW .AND. N.GT.1) W(2) = THEMAX
c set flag if the only error is due to duplicate elements.
      IF (IDUP.EQ.1 .AND. IFLAG.GE.0) IFLAG = -14
      GO TO 70
   60 IF (LP.NE.0) WRITE (LP,99994)
   70 RETURN
99999 FORMAT (39H ERROR RETURN FROM MA28B/BD WITH IFLAG=, I4/I7, 4H ENT,
     * 39HRIES DROPPED FROM STRUCTURE BY MA28A/AD)
99998 FORMAT (36X, 17HN OUT OF RANGE = , I10)
99997 FORMAT (36X, 18HNZ NON POSITIVE = , I10)
99996 FORMAT (36X, 17HLICN TOO SMALL = , I10)
99995 FORMAT (36X, 26HERROR RETURN FROM MA30B/BD)
99994 FORMAT (36H+ERROR RETURN FROM MA28B/BD BECAUSE )
      END
c
c============================================================
c
      SUBROUTINE MA28CD(N, A, LICN, ICN, IKEEP, RHS, W, MTYPE)
c
c this subroutine uses the factors from ma28a/ad or ma28b/bd to
c     solve a system of equations without iterative refinement.
c the parameters are ...
c n   integer  order of matrix  not altered by subroutine.
c a   real/double precision array  length licn.  the same array as
c     was used in the most recent call to ma28a/ad or ma28b/bd.
c licn  integer  length of arrays a and icn.  not altered by
c     subroutine.
c icn    integer array of length licn.  same array as output from
c     ma28a/ad.  unchanged by ma28c/cd.
c ikeep  integer array of length 5*n.  same array as output from
c     ma28a/ad.  unchanged by ma28c/cd.
c rhs    real/double precision array  length n.  on entry, it holds the
c     right hand side.  on exit, the solution vector.
c w      real/double precision array  length n. used as workspace by
c     ma30c/cd.
c mtype  integer  used to tell ma30c/cd to solve the direct equation
c     (mtype=1) or its transpose (mtype.ne.1).
c
      DOUBLE PRECISION A(LICN), RHS(N), W(N), RESID, MRESID, EPS, RMIN
      INTEGER IDISP(2)
      INTEGER ICN(LICN), IKEEP(N,5)
      LOGICAL ABORT1, ABORT2
c common block variables.
c unless otherwise stated common block variables are as in ma28a/ad.
c     those variables referenced by ma28c/cd are mentioned below.
c resid  real/double precision  variable returns maximum residual of
c     equations where pivot was zero.
c mresid  real/double precision variable used by ma28c/cd to
c     communicate between ma28f/fd and ma30h/hd.
c idisp  integer array  length 2  the same as that used by ma28a/ad.
c     it is unchanged by ma28b/bd.
c
c further information on common block variables can be found in block
c     data or ma28a/ad.
      COMMON /MA28FD/ EPS, RMIN, RESID, IRNCP, ICNCP, MINIRN, MINICN,
     * IRANK, ABORT1, ABORT2
      COMMON /MA28GD/ IDISP
      COMMON /MA30HD/ MRESID
c
c this call performs the solution of the set of equations.
      CALL MA30CD(N, ICN, A, LICN, IKEEP, IKEEP(1,4), IKEEP(1,5),
     * IDISP, IKEEP(1,2), IKEEP(1,3), RHS, W, MTYPE)
c transfer common block information.
      RESID = MRESID
      RETURN
      END
      SUBROUTINE MA28DD(N, A, LICN, IVECT, JVECT, NZ, ICN, LENR, LENRL,
     * LENOFF, IP, IQ, IW1, IW, W1, IFLAG)
c this subroutine need never be called by the user directly.
c     it sorts the user's matrix into the structure of the decomposed
c     form and checks for the presence of duplicate entries or
c     non-zeros lying outside the sparsity pattern of the decomposition
c     it also calculates the largest element in the input matrix.
      DOUBLE PRECISION A(LICN), ZERO, W1, AA
      INTEGER IW(N,2), IDISP(2)
      INTEGER ICN(LICN), IVECT(NZ), JVECT(NZ), IP(N), IQ(N),
     * LENR(N), IW1(N,3), LENRL(N), LENOFF(N)
      LOGICAL LBLOCK, GROW, BLOCKL
      COMMON /MA28ED/ LP, MP, LBLOCK, GROW
      COMMON /MA28GD/ IDISP
      DATA ZERO /0.0D0/
      BLOCKL = LENOFF(1).GE.0
c iw1(i,3)  is set to the block in which row i lies and the
c     inverse permutations to ip and iq are set in iw1(.,1) and
c     iw1(.,2) resp.
c pointers to beginning of the part of row i in diagonal and
c   off-diagonal blocks are set in iw(i,2) and iw(i,1) resp.
      IBLOCK = 1
      IW(1,1) = 1
      IW(1,2) = IDISP(1)
      DO 10 I=1,N
        IW1(I,3) = IBLOCK
        IF (IP(I).LT.0) IBLOCK = IBLOCK + 1
        II = IABS(IP(I)+0)
        IW1(II,1) = I
        JJ = IQ(I)
        JJ = IABS(JJ)
        IW1(JJ,2) = I
        IF (I.EQ.1) GO TO 10
        IF (BLOCKL) IW(I,1) = IW(I-1,1) + LENOFF(I-1)
        IW(I,2) = IW(I-1,2) + LENR(I-1)
   10 CONTINUE
c place each non-zero in turn into its correct location
c    in the a/icn array.
      IDISP2 = IDISP(2)
      DO 170 I=1,NZ
c necessary to avoid reference to unassigned element of icn.
        IF (I.GT.IDISP2) GO TO 20
        IF (ICN(I).LT.0) GO TO 170
   20   IOLD = IVECT(I)
        JOLD = JVECT(I)
        AA = A(I)
c this is a dummy loop for following a chain of interchanges.
c   it will be executed nz times in total.
        DO 140 IDUMMY=1,NZ
c perform some validity checks on iold and jold.
          IF (IOLD.LE.N .AND. IOLD.GT.0 .AND. JOLD.LE.N .AND.
     *     JOLD.GT.0) GO TO 30
          IF (LP.NE.0) WRITE (LP,99999) I, A(I), IOLD, JOLD
          IFLAG = -12
          GO TO 180
   30     INEW = IW1(IOLD,1)
          JNEW = IW1(JOLD,2)
c are we in a valid block and is it diagonal or off-diagonal?
          IF (IW1(INEW,3)-IW1(JNEW,3)) 40, 60, 50
   40     IFLAG = -13
          IF (LP.NE.0) WRITE (LP,99998) IOLD, JOLD
          GO TO 180
   50     J1 = IW(INEW,1)
          J2 = J1 + LENOFF(INEW) - 1
          GO TO 110
c element is in diagonal block.
   60     J1 = IW(INEW,2)
          IF (INEW.GT.JNEW) GO TO 70
          J2 = J1 + LENR(INEW) - 1
          J1 = J1 + LENRL(INEW)
          GO TO 110
   70     J2 = J1 + LENRL(INEW)
c binary search of ordered list  .. element in l part of row.
          DO 100 JDUMMY=1,N
            MIDPT = (J1+J2)/2
            JCOMP = IABS(ICN(MIDPT)+0)
            IF (JNEW-JCOMP) 80, 130, 90
   80       J2 = MIDPT
            GO TO 100
   90       J1 = MIDPT
  100     CONTINUE
          IFLAG = -13
          IF (LP.NE.0) WRITE (LP,99997) IOLD, JOLD
          GO TO 180
c linear search ... element in l part of row or off-diagonal blocks.
  110     DO 120 MIDPT=J1,J2
            IF (IABS(ICN(MIDPT)+0).EQ.JNEW) GO TO 130
  120     CONTINUE
          IFLAG = -13
          IF (LP.NE.0) WRITE (LP,99997) IOLD, JOLD
          GO TO 180
c equivalent element of icn is in position midpt.
  130     IF (ICN(MIDPT).LT.0) GO TO 160
          IF (MIDPT.GT.NZ .OR. MIDPT.LE.I) GO TO 150
          W1 = A(MIDPT)
          A(MIDPT) = AA
          AA = W1
          IOLD = IVECT(MIDPT)
          JOLD = JVECT(MIDPT)
          ICN(MIDPT) = -ICN(MIDPT)
  140   CONTINUE
  150   A(MIDPT) = AA
        ICN(MIDPT) = -ICN(MIDPT)
        GO TO 170
  160   A(MIDPT) = A(MIDPT) + AA
c set flag for duplicate elements.
        IFLAG = N + 1
  170 CONTINUE
c reset icn array  and zero elements in l/u but not in a.
c also calculate maximum element of a.
  180 W1 = ZERO
      DO 200 I=1,IDISP2
        IF (ICN(I).LT.0) GO TO 190
        A(I) = ZERO
        GO TO 200
  190   ICN(I) = -ICN(I)
        W1 = DMAX1(W1,DABS(A(I)))
  200 CONTINUE
      RETURN
99999 FORMAT (9H ELEMENT , I6, 12H WITH VALUE , 1PD22.14, 10H HAS INDIC,
     * 3HES , I8, 2H ,, I8/36X, 20HINDICES OUT OF RANGE)
99998 FORMAT (36X, 8HNON-ZERO, I7, 2H ,, I6, 23H IN ZERO OFF-DIAGONAL B,
     * 4HLOCK)
99997 FORMAT (36X, 8H ELEMENT, I6, 2H ,, I6, 23H WAS NOT IN L/U PATTERN)
      END
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias ma30ad
c
c===================================================================
c
      SUBROUTINE MA30AD(NN, ICN, A, LICN, LENR, LENRL, IDISP, IP, IQ,
     * IRN, LIRN, LENC, IFIRST, LASTR, NEXTR, LASTC, NEXTC, IPTR, IPC,
     * U, IFLAG)
c if  the user requires a more convenient data interface then the ma28
c     package should be used.  the ma28 subroutines call the ma30
c     subroutines after checking the user's input data and optionally
c     using mc23a/ad to permute the matrix to block triangular form.
c this package of subroutines (ma30a/ad, ma30b/bd, ma30c/cd and
c     ma30d/dd) performs operations pertinent to the solution of a
c     general sparse n by n system of linear equations (i.e. solve
c     ax=b). structually singular matrices are permitted including
c     those with row or columns consisting entirely of zeros (i.e.
c     including rectangular matrices).  it is assumed that the
c     non-zeros of the matrix a do not differ widely in size.  if
c     necessary a prior call of the scaling subroutine mc19a/ad may be
c     made.
c a discussion of the design of these subroutines is given by duff and
c     reid (acm trans math software 5 pp 18-35,1979 (css 48)) while
c     fuller details of the implementation are given in duff (harwell
c     report aere-r 8730,1977).  the additional pivoting option in
c     ma30a/ad and the use of drop tolerances (see common block
c     ma30i/id) were added to the package after joint work with reid,
c     schaumburg, wasniewski and zlatev (duff, reid, schaumburg,
c     wasniewski and zlatev, harwell report css 135, 1983).
c
c ma30a/ad performs the lu decomposition of the diagonal blocks of the
c     permutation paq of a sparse matrix a, where input permutations
c     p1 and q1 are used to define the diagonal blocks.  there may be
c     non-zeros in the off-diagonal blocks but they are unaffected by
c     ma30a/ad. p and p1 differ only within blocks as do q and q1. the
c     permutations p1 and q1 may be found by calling mc23a/ad or the
c     matrix may be treated as a single block by using p1=q1=i. the
c     matrix non-zeros should be held compactly by rows, although it
c     should be noted that the user can supply the matrix by columns
c     to get the lu decomposition of a transpose.
c
c the parameters are...
c this description should also be consulted for further information on
c     most of the parameters of ma30b/bd and ma30c/cd.
c
c n  is an integer variable which must be set by the user to the order
c     of the matrix.  it is not altered by ma30a/ad.
c icn is an integer array of length licn. positions idisp(2) to
c     licn must be set by the user to contain the column indices of
c     the non-zeros in the diagonal blocks of p1*a*q1. those belonging
c     to a single row must be contiguous but the ordering of column
c     indices with each row is unimportant. the non-zeros of row i
c     precede those of row i+1,i=1,...,n-1 and no wasted space is
c     allowed between the rows.  on output the column indices of the
c     lu decomposition of paq are held in positions idisp(1) to
c     idisp(2), the rows are in pivotal order, and the column indices
c     of the l part of each row are in pivotal order and precede those
c     of u. again there is no wasted space either within a row or
c     between the rows. icn(1) to icn(idisp(1)-1), are neither
c     required nor altered. if mc23a/ad been called, these will hold
c     information about the off-diagonal blocks.
c a is a real/double precision array of length licn whose entries
c     idisp(2) to licn must be set by the user to the  values of the
c     non-zero entries of the matrix in the order indicated by  icn.
c     on output a will hold the lu factors of the matrix where again
c     the position in the matrix is determined by the corresponding
c     values in icn. a(1) to a(idisp(1)-1) are neither required nor
c     altered.
c licn  is an integer variable which must be set by the user to the
c     length of arrays icn and a. it must be big enough for a and icn
c     to hold all the non-zeros of l and u and leave some "elbow
c     room".  it is possible to calculate a minimum value for licn by
c     a preliminary run of ma30a/ad. the adequacy of the elbow room
c     can be judged by the size of the common block variable icncp. it
c     is not altered by ma30a/ad.
c lenr  is an integer array of length n.  on input, lenr(i) should
c     equal the number of non-zeros in row i, i=1,...,n of the
c     diagonal blocks of p1*a*q1. on output, lenr(i) will equal the
c     total number of non-zeros in row i of l and row i of u.
c lenrl  is an integer array of length n. on output from ma30a/ad,
c     lenrl(i) will hold the number of non-zeros in row i of l.
c idisp  is an integer array of length 2. the user should set idisp(1)
c     to be the first available position in a/icn for the lu
c     decomposition while idisp(2) is set to the position in a/icn of
c     the first non-zero in the diagonal blocks of p1*a*q1. on output,
c     idisp(1) will be unaltered while idisp(2) will be set to the
c     position in a/icn of the last non-zero of the lu decomposition.
c ip  is an integer array of length n which holds a permutation of
c     the integers 1 to n.  on input to ma30a/ad, the absolute value of
c     ip(i) must be set to the row of a which is row i of p1*a*q1. a
c     negative value for ip(i) indicates that row i is at the end of a
c     diagonal block.  on output from ma30a/ad, ip(i) indicates the row
c     of a which is the i th row in paq. ip(i) will still be negative
c     for the last row of each block (except the last).
c iq is an integer array of length n which again holds a
c     permutation of the integers 1 to n.  on input to ma30a/ad, iq(j)
c     must be set to the column of a which is column j of p1*a*q1. on
c     output from ma30a/ad, the absolute value of iq(j) indicates the
c     column of a which is the j th in paq.  for rows, i say, in which
c     structural or numerical singularity is detected iq(i) is
c     negated.
c irn  is an integer array of length lirn used as workspace by
c     ma30a/ad.
c lirn  is an integer variable. it should be greater than the
c     largest number of non-zeros in a diagonal block of p1*a*q1 but
c     need not be as large as licn. it is the length of array irn and
c     should be large enough to hold the active part of any block,
c     plus some "elbow room", the  a posteriori  adequacy of which can
c     be estimated by examining the size of common block variable
c     irncp.
c lenc,ifirst,lastr,nextr,lastc,nextc are all integer arrays of
c     length n which are used as workspace by ma30a/ad.  if nsrch is
c     set to a value less than or equal to n, then arrays lastc and
c     nextc are not referenced by ma30a/ad and so can be dummied in
c     the call to ma30a/ad.
c iptr,ipc are integer arrays of length n which are used as workspace
c     by ma30a/ad.
c u  is a real/double precision variable which should be set by the
c     user to a value between 0. and 1.0. if less than zero it is
c     reset to zero and if its value is 1.0 or greater it is reset to
c     0.9999 (0.999999999 in d version).  it determines the balance
c     between pivoting for sparsity and for stability, values near
c     zero emphasizing sparsity and values near one emphasizing
c     stability. we recommend u=0.1 as a posible first trial value.
c     the stability can be judged by a later call to mc24a/ad or by
c     setting lbig to .true.
c iflag  is an integer variable. it will have a non-negative value if
c     ma30a/ad is successful. negative values indicate error
c     conditions while positive values indicate that the matrix has
c     been successfully decomposed but is singular. for each non-zero
c     value, an appropriate message is output on unit lp.  possible
c     non-zero values for iflag are ...
c
c -1  the matrix is structually singular with rank given by irank in
c     common block ma30f/fd.
c +1  if, however, the user wants the lu decomposition of a
c     structurally singular matrix and sets the common block variable
c     abort1 to .false., then, in the event of singularity and a
c     successful decomposition, iflag is returned with the value +1
c     and no message is output.
c -2  the matrix is numerically singular (it may also be structually
c     singular) with estimated rank given by irank in common block
c     ma30f/fd.
c +2  the  user can choose to continue the decomposition even when a
c     zero pivot is encountered by setting common block variable
c     abort2 to .false.  if a singularity is encountered, iflag will
c     then return with a value of +2, and no message is output if the
c     decomposition has been completed successfully.
c -3  lirn has not been large enough to continue with the
c     decomposition.  if the stage was zero then common block variable
c     minirn gives the length sufficient to start the decomposition on
c     this block.  for a successful decomposition on this block the
c     user should make lirn slightly (say about n/2) greater than this
c     value.
c -4  licn not large enough to continue with the decomposition.
c -5  the decomposition has been completed but some of the lu factors
c     have been discarded to create enough room in a/icn to continue
c     the decomposition. the variable minicn in common block ma30f/fd
c     then gives the size that licn should be to enable the
c     factorization to be successful.  if the user sets common block
c     variable abort3 to .true., then the subroutine will exit
c     immediately instead of destroying any factors and continuing.
c -6  both licn and lirn are too small. termination has been caused by
c     lack of space in irn (see error iflag= -3), but already some of
c     the lu factors in a/icn have been lost (see error iflag= -5).
c     minicn gives the minimum amount of space required in a/icn for
c     decomposition up to this point.
c
      DOUBLE PRECISION A(LICN), U, AU, UMAX, AMAX, ZERO, PIVRAT, PIVR,
     * TOL, BIG, ANEW, AANEW, SCALE
      INTEGER IPTR(NN), PIVOT, PIVEND, DISPC, OLDPIV, OLDEND, PIVROW,
     * ROWI, IPC(NN), IDISP(2), COLUPD
      INTEGER ICN(LICN), LENR(NN), LENRL(NN), IP(NN), IQ(NN),
     * LENC(NN), IRN(LIRN), IFIRST(NN), LASTR(NN), NEXTR(NN),
     * LASTC(NN), NEXTC(NN)
      LOGICAL ABORT1, ABORT2, ABORT3, LBIG
c for comments of common block variables see block data subprogram.
      COMMON /MA30ED/ LP, ABORT1, ABORT2, ABORT3
      COMMON /MA30FD/ IRNCP, ICNCP, IRANK, MINIRN, MINICN
      COMMON /MA30ID/ TOL, BIG, NDROP, NSRCH, LBIG
      COMMON /LPIVOT/ LPIV(10),LNPIV(10),MAPIV,MANPIV,IAVPIV,
     *                IANPIV,KOUNTL
c
      DATA UMAX/.999999999D0/
      DATA ZERO /0.0D0/
      MSRCH = NSRCH
      NDROP = 0
      DO 1272 KK=1,10
        LNPIV(KK)=0
        LPIV(KK)=0
 1272 CONTINUE
      MAPIV = 0
      MANPIV = 0
      IAVPIV = 0
      IANPIV = 0
      KOUNTL = 0
      MINIRN = 0
      MINICN = IDISP(1) - 1
      MOREI = 0
      IRANK = NN
      IRNCP = 0
      ICNCP = 0
      IFLAG = 0
c reset u if necessary.
      U = DMIN1(U,UMAX)
c ibeg is the position of the next pivot row after elimination step
c     using it.
      U = DMAX1(U,ZERO)
      IBEG = IDISP(1)
c iactiv is the position of the first entry in the active part of a/icn.
      IACTIV = IDISP(2)
c nzrow is current number of non-zeros in active and unprocessed part
c     of row file icn.
      NZROW = LICN - IACTIV + 1
      MINICN = NZROW + MINICN
c
c count the number of diagonal blocks and set up pointers to the
c     beginnings of the rows.
c num is the number of diagonal blocks.
      NUM = 1
      IPTR(1) = IACTIV
      IF (NN.EQ.1) GO TO 20
      NNM1 = NN - 1
      DO 10 I=1,NNM1
        IF (IP(I).LT.0) NUM = NUM + 1
        IPTR(I+1) = IPTR(I) + LENR(I)
   10 CONTINUE
c ilast is the last row in the previous block.
   20 ILAST = 0
c
c ***********************************************
c ****    lu decomposition of block nblock   ****
c ***********************************************
c
c each pass through this loop performs lu decomposition on one
c     of the diagonal blocks.
      DO 1000 NBLOCK=1,NUM
        ISTART = ILAST + 1
        DO 30 IROWS=ISTART,NN
          IF (IP(IROWS).LT.0) GO TO 40
   30   CONTINUE
        IROWS = NN
   40   ILAST = IROWS
c n is the number of rows in the current block.
c istart is the index of the first row in the current block.
c ilast is the index of the last row in the current block.
c iactiv is the position of the first entry in the block.
c itop is the position of the last entry in the block.
        N = ILAST - ISTART + 1
        IF (N.NE.1) GO TO 90
c
c code for dealing with 1x1 block.
        LENRL(ILAST) = 0
        ISING = ISTART
        IF (LENR(ILAST).NE.0) GO TO 50
c block is structurally singular.
        IRANK = IRANK - 1
        ISING = -ISING
        IF (IFLAG.NE.2 .AND. IFLAG.NE.-5) IFLAG = 1
        IF (.NOT.ABORT1) GO TO 80
        IDISP(2) = IACTIV
        IFLAG = -1
        IF (LP.NE.0) WRITE (LP,99999)
c     return
        GO TO 1120
   50   SCALE = DABS(A(IACTIV))
        IF (SCALE.EQ.ZERO) GO TO 60
        IF (LBIG) BIG = DMAX1(BIG,SCALE)
        GO TO 70
   60   ISING = -ISING
        IRANK = IRANK - 1
        IPTR(ILAST) = 0
        IF (IFLAG.NE.-5) IFLAG = 2
        IF (.NOT.ABORT2) GO TO 70
        IDISP(2) = IACTIV
        IFLAG = -2
        IF (LP.NE.0) WRITE (LP,99998)
        GO TO 1120
   70   A(IBEG) = A(IACTIV)
        ICN(IBEG) = ICN(IACTIV)
        IACTIV = IACTIV + 1
        IPTR(ISTART) = 0
        IBEG = IBEG + 1
        NZROW = NZROW - 1
   80   LASTR(ISTART) = ISTART
        IPC(ISTART) = -ISING
        GO TO 1000
c
c non-trivial block.
   90   ITOP = LICN
        IF (ILAST.NE.NN) ITOP = IPTR(ILAST+1) - 1
c
c set up column oriented storage.
        DO 100 I=ISTART,ILAST
          LENRL(I) = 0
          LENC(I) = 0
  100   CONTINUE
        IF (ITOP-IACTIV.LT.LIRN) GO TO 110
        MINIRN = ITOP - IACTIV + 1
        PIVOT = ISTART - 1
        GO TO 1100
c
c calculate column counts.
  110   DO 120 II=IACTIV,ITOP
          I = ICN(II)
          LENC(I) = LENC(I) + 1
  120   CONTINUE
c set up column pointers so that ipc(j) points to position after end
c     of column j in column file.
        IPC(ILAST) = LIRN + 1
        J1 = ISTART + 1
        DO 130 JJ=J1,ILAST
          J = ILAST - JJ + J1 - 1
          IPC(J) = IPC(J+1) - LENC(J+1)
  130   CONTINUE
        DO 150 INDROW=ISTART,ILAST
          J1 = IPTR(INDROW)
          J2 = J1 + LENR(INDROW) - 1
          IF (J1.GT.J2) GO TO 150
          DO 140 JJ=J1,J2
            J = ICN(JJ)
            IPOS = IPC(J) - 1
            IRN(IPOS) = INDROW
            IPC(J) = IPOS
  140     CONTINUE
  150   CONTINUE
c dispc is the lowest indexed active location in the column file.
        DISPC = IPC(ISTART)
        NZCOL = LIRN - DISPC + 1
        MINIRN = MAX0(NZCOL,MINIRN)
        NZMIN = 1
c
c initialize array ifirst.  ifirst(i) = +/- k indicates that row/col
c     k has i non-zeros.  if ifirst(i) = 0, there is no row or column
c     with i non zeros.
        DO 160 I=1,N
          IFIRST(I) = 0
  160   CONTINUE
c
c compute ordering of row and column counts.
c first run through columns (from column n to column 1).
        DO 180 JJ=ISTART,ILAST
          J = ILAST - JJ + ISTART
          NZ = LENC(J)
          IF (NZ.NE.0) GO TO 170
          IPC(J) = 0
          GO TO 180
  170     IF (NSRCH.LE.NN) GO TO 180
          ISW = IFIRST(NZ)
          IFIRST(NZ) = -J
          LASTC(J) = 0
          NEXTC(J) = -ISW
          ISW1 = IABS(ISW)
          IF (ISW.NE.0) LASTC(ISW1) = J
  180   CONTINUE
c now run through rows (again from n to 1).
        DO 210 II=ISTART,ILAST
          I = ILAST - II + ISTART
          NZ = LENR(I)
          IF (NZ.NE.0) GO TO 190
          IPTR(I) = 0
          LASTR(I) = 0
          GO TO 210
  190     ISW = IFIRST(NZ)
          IFIRST(NZ) = I
          IF (ISW.GT.0) GO TO 200
          NEXTR(I) = 0
          LASTR(I) = ISW
          GO TO 210
  200     NEXTR(I) = ISW
          LASTR(I) = LASTR(ISW)
          LASTR(ISW) = I
  210   CONTINUE
c
c
c **********************************************
c ****    start of main elimination loop    ****
c **********************************************
        DO 980 PIVOT=ISTART,ILAST
c
c first find the pivot using markowitz criterion with stability
c     control.
c jcost is the markowitz cost of the best pivot so far,.. this
c     pivot is in row ipiv and column jpiv.
          NZ2 = NZMIN
          JCOST = N*N
c
c examine rows/columns in order of ascending count.
          DO 340 L=1,2
            PIVRAT = ZERO
            ISRCH = 1
            LL = L
c a pass with l equal to 2 is only performed in the case of singularity.
            DO 330 NZ=NZ2,N
              IF (JCOST.LE.(NZ-1)**2) GO TO 420
              IJFIR = IFIRST(NZ)
              IF (IJFIR) 230, 220, 240
  220         IF (LL.EQ.1) NZMIN = NZ + 1
              GO TO 330
  230         LL = 2
              IJFIR = -IJFIR
              GO TO 290
  240         LL = 2
c scan rows with nz non-zeros.
              DO 270 IDUMMY=1,N
                IF (JCOST.LE.(NZ-1)**2) GO TO 420
                IF (ISRCH.GT.MSRCH) GO TO 420
                IF (IJFIR.EQ.0) GO TO 280
c row ijfir is now examined.
                I = IJFIR
                IJFIR = NEXTR(I)
c first calculate multiplier threshold level.
                AMAX = ZERO
                J1 = IPTR(I) + LENRL(I)
                J2 = IPTR(I) + LENR(I) - 1
                DO 250 JJ=J1,J2
                  AMAX = DMAX1(AMAX,DABS(A(JJ)))
  250           CONTINUE
                AU = AMAX*U
                ISRCH = ISRCH + 1
c scan row for possible pivots
                DO 260 JJ=J1,J2
                  IF (DABS(A(JJ)).LE.AU .AND. L.EQ.1) GO TO 260
                  J = ICN(JJ)
                  KCOST = (NZ-1)*(LENC(J)-1)
                  IF (KCOST.GT.JCOST) GO TO 260
                  PIVR = ZERO
                  IF (AMAX.NE.ZERO) PIVR = DABS(A(JJ))/AMAX
                  IF (KCOST.EQ.JCOST .AND. (PIVR.LE.PIVRAT .OR.
     *             NSRCH.GT.NN+1)) GO TO 260
c best pivot so far is found.
                  JCOST = KCOST
                  IJPOS = JJ
                  IPIV = I
                  JPIV = J
                  IF (MSRCH.GT.NN+1 .AND. JCOST.LE.(NZ-1)**2) GO TO 420
                  PIVRAT = PIVR
  260           CONTINUE
  270         CONTINUE
c
c columns with nz non-zeros now examined.
  280         IJFIR = IFIRST(NZ)
              IJFIR = -LASTR(IJFIR)
  290         IF (JCOST.LE.NZ*(NZ-1)) GO TO 420
              IF (MSRCH.LE.NN) GO TO 330
              DO 320 IDUMMY=1,N
                IF (IJFIR.EQ.0) GO TO 330
                J = IJFIR
                IJFIR = NEXTC(IJFIR)
                I1 = IPC(J)
                I2 = I1 + NZ - 1
c scan column j.
                DO 310 II=I1,I2
                  I = IRN(II)
                  KCOST = (NZ-1)*(LENR(I)-LENRL(I)-1)
                  IF (KCOST.GE.JCOST) GO TO 310
c pivot has best markowitz count so far ... now check its
c     suitability on numeric grounds by examining the other non-zeros
c     in its row.
                  J1 = IPTR(I) + LENRL(I)
                  J2 = IPTR(I) + LENR(I) - 1
c we need a stability check on singleton columns because of possible
c     problems with underdetermined systems.
                  AMAX = ZERO
                  DO 300 JJ=J1,J2
                    AMAX = DMAX1(AMAX,DABS(A(JJ)))
                    IF (ICN(JJ).EQ.J) JPOS = JJ
  300             CONTINUE
                  IF (DABS(A(JPOS)).LE.AMAX*U .AND. L.EQ.1) GO TO 310
                  JCOST = KCOST
                  IPIV = I
                  JPIV = J
                  IJPOS = JPOS
                  IF (AMAX.NE.ZERO) PIVRAT = DABS(A(JPOS))/AMAX
                  IF (JCOST.LE.NZ*(NZ-1)) GO TO 420
  310           CONTINUE
c
  320         CONTINUE
c
  330       CONTINUE
c in the event of singularity, we must make sure all rows and columns
c are tested.
            MSRCH = N
c
c matrix is numerically or structurally singular  ... which it is will
c     be diagnosed later.
            IRANK = IRANK - 1
  340     CONTINUE
c assign rest of rows and columns to ordering array.
c matrix is structurally singular.
          IF (IFLAG.NE.2 .AND. IFLAG.NE.-5) IFLAG = 1
          IRANK = IRANK - ILAST + PIVOT + 1
          IF (.NOT.ABORT1) GO TO 350
          IDISP(2) = IACTIV
          IFLAG = -1
          IF (LP.NE.0) WRITE (LP,99999)
          GO TO 1120
  350     K = PIVOT - 1
          DO 390 I=ISTART,ILAST
            IF (LASTR(I).NE.0) GO TO 390
            K = K + 1
            LASTR(I) = K
            IF (LENRL(I).EQ.0) GO TO 380
            MINICN = MAX0(MINICN,NZROW+IBEG-1+MOREI+LENRL(I))
            IF (IACTIV-IBEG.GE.LENRL(I)) GO TO 360
            CALL MA30DD(A, ICN, IPTR(ISTART), N, IACTIV, ITOP, .TRUE.)
c check now to see if ma30d/dd has created enough available space.
            IF (IACTIV-IBEG.GE.LENRL(I)) GO TO 360
c create more space by destroying previously created lu factors.
            MOREI = MOREI + IBEG - IDISP(1)
            IBEG = IDISP(1)
            IF (LP.NE.0) WRITE (LP,99997)
            IFLAG = -5
            IF (ABORT3) GO TO 1090
  360       J1 = IPTR(I)
            J2 = J1 + LENRL(I) - 1
            IPTR(I) = 0
            DO 370 JJ=J1,J2
              A(IBEG) = A(JJ)
              ICN(IBEG) = ICN(JJ)
              ICN(JJ) = 0
              IBEG = IBEG + 1
  370       CONTINUE
            NZROW = NZROW - LENRL(I)
  380       IF (K.EQ.ILAST) GO TO 400
  390     CONTINUE
  400     K = PIVOT - 1
          DO 410 I=ISTART,ILAST
            IF (IPC(I).NE.0) GO TO 410
            K = K + 1
            IPC(I) = K
            IF (K.EQ.ILAST) GO TO 990
  410     CONTINUE
c
c the pivot has now been found in position (ipiv,jpiv) in location
c     ijpos in row file.
c update column and row ordering arrays to correspond with removal
c     of the active part of the matrix.
  420     ISING = PIVOT
          IF (A(IJPOS).NE.ZERO) GO TO 430
c numerical singularity is recorded here.
          ISING = -ISING
          IF (IFLAG.NE.-5) IFLAG = 2
          IF (.NOT.ABORT2) GO TO 430
          IDISP(2) = IACTIV
          IFLAG = -2
          IF (LP.NE.0) WRITE (LP,99998)
          GO TO 1120
  430     OLDPIV = IPTR(IPIV) + LENRL(IPIV)
          OLDEND = IPTR(IPIV) + LENR(IPIV) - 1
c changes to column ordering.
          IF (NSRCH.LE.NN) GO TO 460
          COLUPD = NN + 1
            LENPP = OLDEND-OLDPIV+1
            IF (LENPP.LT.4) LPIV(1) = LPIV(1) + 1
            IF (LENPP.GE.4 .AND. LENPP.LE.6) LPIV(2) = LPIV(2) + 1
            IF (LENPP.GE.7 .AND. LENPP.LE.10) LPIV(3) = LPIV(3) + 1
            IF (LENPP.GE.11 .AND. LENPP.LE.15) LPIV(4) = LPIV(4) + 1
            IF (LENPP.GE.16 .AND. LENPP.LE.20) LPIV(5) = LPIV(5) + 1
            IF (LENPP.GE.21 .AND. LENPP.LE.30) LPIV(6) = LPIV(6) + 1
            IF (LENPP.GE.31 .AND. LENPP.LE.50) LPIV(7) = LPIV(7) + 1
            IF (LENPP.GE.51 .AND. LENPP.LE.70) LPIV(8) = LPIV(8) + 1
            IF (LENPP.GE.71 .AND. LENPP.LE.100) LPIV(9) = LPIV(9) + 1
            IF (LENPP.GE.101) LPIV(10) = LPIV(10) + 1
            MAPIV = MAX0(MAPIV,LENPP)
            IAVPIV = IAVPIV + LENPP
          DO 450 JJ=OLDPIV,OLDEND
            J = ICN(JJ)
            LC = LASTC(J)
            NC = NEXTC(J)
            NEXTC(J) = -COLUPD
            IF (JJ.NE.IJPOS) COLUPD = J
            IF (NC.NE.0) LASTC(NC) = LC
            IF (LC.EQ.0) GO TO 440
            NEXTC(LC) = NC
            GO TO 450
  440       NZ = LENC(J)
            ISW = IFIRST(NZ)
            IF (ISW.GT.0) LASTR(ISW) = -NC
            IF (ISW.LT.0) IFIRST(NZ) = -NC
  450     CONTINUE
c changes to row ordering.
  460     I1 = IPC(JPIV)
          I2 = I1 + LENC(JPIV) - 1
          DO 480 II=I1,I2
            I = IRN(II)
            LR = LASTR(I)
            NR = NEXTR(I)
            IF (NR.NE.0) LASTR(NR) = LR
            IF (LR.LE.0) GO TO 470
            NEXTR(LR) = NR
            GO TO 480
  470       NZ = LENR(I) - LENRL(I)
            IF (NR.NE.0) IFIRST(NZ) = NR
            IF (NR.EQ.0) IFIRST(NZ) = LR
  480     CONTINUE
c
c move pivot to position lenrl+1 in pivot row and move pivot row
c     to the beginning of the available storage.
c the l part and the pivot in the old copy of the pivot row is
c     nullified while, in the strictly upper triangular part, the
c     column indices, j say, are overwritten by the corresponding
c     entry of iq (iq(j)) and iq(j) is set to the negative of the
c     displacement of the column index from the pivot entry.
          IF (OLDPIV.EQ.IJPOS) GO TO 490
          AU = A(OLDPIV)
          A(OLDPIV) = A(IJPOS)
          A(IJPOS) = AU
          ICN(IJPOS) = ICN(OLDPIV)
          ICN(OLDPIV) = JPIV
c check to see if there is space immediately available in a/icn to
c     hold new copy of pivot row.
  490     MINICN = MAX0(MINICN,NZROW+IBEG-1+MOREI+LENR(IPIV))
          IF (IACTIV-IBEG.GE.LENR(IPIV)) GO TO 500
          CALL MA30DD(A, ICN, IPTR(ISTART), N, IACTIV, ITOP, .TRUE.)
          OLDPIV = IPTR(IPIV) + LENRL(IPIV)
          OLDEND = IPTR(IPIV) + LENR(IPIV) - 1
c check now to see if ma30d/dd has created enough available space.
          IF (IACTIV-IBEG.GE.LENR(IPIV)) GO TO 500
c create more space by destroying previously created lu factors.
          MOREI = MOREI + IBEG - IDISP(1)
          IBEG = IDISP(1)
          IF (LP.NE.0) WRITE (LP,99997)
          IFLAG = -5
          IF (ABORT3) GO TO 1090
          IF (IACTIV-IBEG.GE.LENR(IPIV)) GO TO 500
c there is still not enough room in a/icn.
          IFLAG = -4
          GO TO 1090
c copy pivot row and set up iq array.
  500     IJPOS = 0
          J1 = IPTR(IPIV)
c
          DO 530 JJ=J1,OLDEND
            A(IBEG) = A(JJ)
            ICN(IBEG) = ICN(JJ)
            IF (IJPOS.NE.0) GO TO 510
            IF (ICN(JJ).EQ.JPIV) IJPOS = IBEG
            ICN(JJ) = 0
            GO TO 520
  510       K = IBEG - IJPOS
            J = ICN(JJ)
            ICN(JJ) = IQ(J)
            IQ(J) = -K
  520       IBEG = IBEG + 1
  530     CONTINUE
c
          IJP1 = IJPOS + 1
          PIVEND = IBEG - 1
          LENPIV = PIVEND - IJPOS
          NZROW = NZROW - LENRL(IPIV) - 1
          IPTR(IPIV) = OLDPIV + 1
          IF (LENPIV.EQ.0) IPTR(IPIV) = 0
c
c remove pivot row (including pivot) from column oriented file.
          DO 560 JJ=IJPOS,PIVEND
            J = ICN(JJ)
            I1 = IPC(J)
            LENC(J) = LENC(J) - 1
c i2 is last position in new column.
            I2 = IPC(J) + LENC(J) - 1
            IF (I2.LT.I1) GO TO 550
            DO 540 II=I1,I2
              IF (IRN(II).NE.IPIV) GO TO 540
              IRN(II) = IRN(I2+1)
              GO TO 550
  540       CONTINUE
  550       IRN(I2+1) = 0
  560     CONTINUE
          NZCOL = NZCOL - LENPIV - 1
c
c go down the pivot column and for each row with a non-zero add
c     the appropriate multiple of the pivot row to it.
c we loop on the number of non-zeros in the pivot column since
c     ma30d/dd may change its actual position.
c
          NZPC = LENC(JPIV)
          IF (NZPC.EQ.0) GO TO 900
          DO 840 III=1,NZPC
            II = IPC(JPIV) + III - 1
            I = IRN(II)
c search row i for non-zero to be eliminated, calculate multiplier,
c     and place it in position lenrl+1 in its row.
c  idrop is the number of non-zero entries dropped from row    i
c        because these fall beneath tolerance level.
c
            IDROP = 0
            J1 = IPTR(I) + LENRL(I)
            IEND = IPTR(I) + LENR(I) - 1
            DO 570 JJ=J1,IEND
              IF (ICN(JJ).NE.JPIV) GO TO 570
c if pivot is zero, rest of column is and so multiplier is zero.
              AU = ZERO
              IF (A(IJPOS).NE.ZERO) AU = -A(JJ)/A(IJPOS)
              IF (LBIG) BIG = DMAX1(BIG,DABS(AU))
              A(JJ) = A(J1)
              A(J1) = AU
              ICN(JJ) = ICN(J1)
              ICN(J1) = JPIV
              LENRL(I) = LENRL(I) + 1
              GO TO 580
  570       CONTINUE
c jump if pivot row is a singleton.
  580       IF (LENPIV.EQ.0) GO TO 840
c now perform necessary operations on rest of non-pivot row i.
            ROWI = J1 + 1
            IOP = 0
c jump if all the pivot row causes fill-in.
            IF (ROWI.GT.IEND) GO TO 650
c perform operations on current non-zeros in row i.
c innermost loop.
            LENPP = IEND-ROWI+1
            IF (LENPP.LT.4) LNPIV(1) = LNPIV(1) + 1
            IF (LENPP.GE.4 .AND. LENPP.LE.6) LNPIV(2) = LNPIV(2) + 1
            IF (LENPP.GE.7 .AND. LENPP.LE.10) LNPIV(3) = LNPIV(3) + 1
            IF (LENPP.GE.11 .AND. LENPP.LE.15) LNPIV(4) = LNPIV(4) + 1
            IF (LENPP.GE.16 .AND. LENPP.LE.20) LNPIV(5) = LNPIV(5) + 1
            IF (LENPP.GE.21 .AND. LENPP.LE.30) LNPIV(6) = LNPIV(6) + 1
            IF (LENPP.GE.31 .AND. LENPP.LE.50) LNPIV(7) = LNPIV(7) + 1
            IF (LENPP.GE.51 .AND. LENPP.LE.70) LNPIV(8) = LNPIV(8) + 1
            IF (LENPP.GE.71 .AND. LENPP.LE.100) LNPIV(9) = LNPIV(9) + 1
            IF (LENPP.GE.101) LNPIV(10) = LNPIV(10) + 1
            MANPIV = MAX0(MANPIV,LENPP)
            IANPIV = IANPIV + LENPP
            KOUNTL = KOUNTL + 1
            DO 590 JJ=ROWI,IEND
              J = ICN(JJ)
              IF (IQ(J).GT.0) GO TO 590
              IOP = IOP + 1
              PIVROW = IJPOS - IQ(J)
              A(JJ) = A(JJ) + AU*A(PIVROW)
              IF (LBIG) BIG = DMAX1(DABS(A(JJ)),BIG)
              ICN(PIVROW) = -ICN(PIVROW)
              IF (DABS(A(JJ)).LT.TOL) IDROP = IDROP + 1
  590       CONTINUE
c
c  jump if no non-zeros in non-pivot row have been removed
c       because these are beneath the drop-tolerance  tol.
c
            IF (IDROP.EQ.0) GO TO 650
c
c  run through non-pivot row compressing row so that only
c      non-zeros greater than   tol   are stored.  all non-zeros
c      less than   tol   are also removed from the column structure.
c
            JNEW = ROWI
            DO 630 JJ=ROWI,IEND
              IF (DABS(A(JJ)).LT.TOL) GO TO 600
              A(JNEW) = A(JJ)
              ICN(JNEW) = ICN(JJ)
              JNEW = JNEW + 1
              GO TO 630
c
c  remove non-zero entry from column structure.
c
  600         J = ICN(JJ)
              I1 = IPC(J)
              I2 = I1 + LENC(J) - 1
              DO 610 II=I1,I2
                IF (IRN(II).EQ.I) GO TO 620
  610         CONTINUE
  620         IRN(II) = IRN(I2)
              IRN(I2) = 0
              LENC(J) = LENC(J) - 1
              IF (NSRCH.LE.NN) GO TO 630
c remove column from column chain and place in update chain.
              IF (NEXTC(J).LT.0) GO TO 630
c jump if column already in update chain.
              LC = LASTC(J)
              NC = NEXTC(J)
              NEXTC(J) = -COLUPD
              COLUPD = J
              IF (NC.NE.0) LASTC(NC) = LC
              IF (LC.EQ.0) GO TO 622
              NEXTC(LC) = NC
              GO TO 630
  622         NZ = LENC(J) + 1
              ISW = IFIRST(NZ)
              IF (ISW.GT.0) LASTR(ISW) = -NC
              IF (ISW.LT.0) IFIRST(NZ) = -NC
  630       CONTINUE
            DO 640 JJ=JNEW,IEND
              ICN(JJ) = 0
  640       CONTINUE
c the value of idrop might be different from that calculated earlier
c     because, we may now have dropped some non-zeros which were not
c     modified by the pivot row.
            IDROP = IEND + 1 - JNEW
            IEND = JNEW - 1
            LENR(I) = LENR(I) - IDROP
            NZROW = NZROW - IDROP
            NZCOL = NZCOL - IDROP
            NDROP = NDROP + IDROP
  650       IFILL = LENPIV - IOP
c jump is if there is no fill-in.
            IF (IFILL.EQ.0) GO TO 750
c now for the fill-in.
            MINICN = MAX0(MINICN,MOREI+IBEG-1+NZROW+IFILL+LENR(I))
c see if there is room for fill-in.
c get maximum space for row i in situ.
            DO 660 JDIFF=1,IFILL
              JNPOS = IEND + JDIFF
              IF (JNPOS.GT.LICN) GO TO 670
              IF (ICN(JNPOS).NE.0) GO TO 670
  660       CONTINUE
c there is room for all the fill-in after the end of the row so it
c     can be left in situ.
c next available space for fill-in.
            IEND = IEND + 1
            GO TO 750
c jmore spaces for fill-in are required in front of row.
  670       JMORE = IFILL - JDIFF + 1
            I1 = IPTR(I)
c we now look in front of the row to see if there is space for
c     the rest of the fill-in.
            DO 680 JDIFF=1,JMORE
              JNPOS = I1 - JDIFF
              IF (JNPOS.LT.IACTIV) GO TO 690
              IF (ICN(JNPOS).NE.0) GO TO 700
  680       CONTINUE
  690       JNPOS = I1 - JMORE
            GO TO 710
c whole row must be moved to the beginning of available storage.
  700       JNPOS = IACTIV - LENR(I) - IFILL
c jump if there is space immediately available for the shifted row.
  710       IF (JNPOS.GE.IBEG) GO TO 730
            CALL MA30DD(A, ICN, IPTR(ISTART), N, IACTIV, ITOP, .TRUE.)
            I1 = IPTR(I)
            IEND = I1 + LENR(I) - 1
            JNPOS = IACTIV - LENR(I) - IFILL
            IF (JNPOS.GE.IBEG) GO TO 730
c no space available so try to create some by throwing away previous
c     lu decomposition.
            MOREI = MOREI + IBEG - IDISP(1) - LENPIV - 1
            IF (LP.NE.0) WRITE (LP,99997)
            IFLAG = -5
            IF (ABORT3) GO TO 1090
c keep record of current pivot row.
            IBEG = IDISP(1)
            ICN(IBEG) = JPIV
            A(IBEG) = A(IJPOS)
            IJPOS = IBEG
            DO 720 JJ=IJP1,PIVEND
              IBEG = IBEG + 1
              A(IBEG) = A(JJ)
              ICN(IBEG) = ICN(JJ)
  720       CONTINUE
            IJP1 = IJPOS + 1
            PIVEND = IBEG
            IBEG = IBEG + 1
            IF (JNPOS.GE.IBEG) GO TO 730
c this still does not give enough room.
            IFLAG = -4
            GO TO 1090
  730       IACTIV = MIN0(IACTIV,JNPOS)
c move non-pivot row i.
            IPTR(I) = JNPOS
            DO 740 JJ=I1,IEND
              A(JNPOS) = A(JJ)
              ICN(JNPOS) = ICN(JJ)
              JNPOS = JNPOS + 1
              ICN(JJ) = 0
  740       CONTINUE
c first new available space.
            IEND = JNPOS
  750       NZROW = NZROW + IFILL
c innermost fill-in loop which also resets icn.
            IDROP = 0
            DO 830 JJ=IJP1,PIVEND
              J = ICN(JJ)
              IF (J.LT.0) GO TO 820
              ANEW = AU*A(JJ)
              AANEW = DABS(ANEW)
              IF (AANEW.GE.TOL) GO TO 760
              IDROP = IDROP + 1
              NDROP = NDROP + 1
              NZROW = NZROW - 1
              MINICN = MINICN - 1
              IFILL = IFILL - 1
              GO TO 830
  760         IF (LBIG) BIG = DMAX1(AANEW,BIG)
              A(IEND) = ANEW
              ICN(IEND) = J
              IEND = IEND + 1
c
c put new entry in column file.
              MINIRN = MAX0(MINIRN,NZCOL+LENC(J)+1)
              JEND = IPC(J) + LENC(J)
              JROOM = NZPC - III + 1 + LENC(J)
              IF (JEND.GT.LIRN) GO TO 770
              IF (IRN(JEND).EQ.0) GO TO 810
  770         IF (JROOM.LT.DISPC) GO TO 780
c compress column file to obtain space for new copy of column.
              CALL MA30DD(A, IRN, IPC(ISTART), N, DISPC, LIRN, .FALSE.)
              IF (JROOM.LT.DISPC) GO TO 780
              JROOM = DISPC - 1
              IF (JROOM.GE.LENC(J)+1) GO TO 780
c column file is not large enough.
              GO TO 1100
c copy column to beginning of file.
  780         JBEG = IPC(J)
              JEND = IPC(J) + LENC(J) - 1
              JZERO = DISPC - 1
              DISPC = DISPC - JROOM
              IDISPC = DISPC
              DO 790 II=JBEG,JEND
                IRN(IDISPC) = IRN(II)
                IRN(II) = 0
                IDISPC = IDISPC + 1
  790         CONTINUE
              IPC(J) = DISPC
              JEND = IDISPC
              DO 800 II=JEND,JZERO
                IRN(II) = 0
  800         CONTINUE
  810         IRN(JEND) = I
              NZCOL = NZCOL + 1
              LENC(J) = LENC(J) + 1
c end of adjustment to column file.
              GO TO 830
c
  820         ICN(JJ) = -J
  830       CONTINUE
            IF (IDROP.EQ.0) GO TO 834
            DO 832 KDROP=1,IDROP
            ICN(IEND) = 0
            IEND = IEND + 1
  832       CONTINUE
  834       LENR(I) = LENR(I) + IFILL
c end of scan of pivot column.
  840     CONTINUE
c
c
c remove pivot column from column oriented storage and update row
c     ordering arrays.
          I1 = IPC(JPIV)
          I2 = IPC(JPIV) + LENC(JPIV) - 1
          NZCOL = NZCOL - LENC(JPIV)
          DO 890 II=I1,I2
            I = IRN(II)
            IRN(II) = 0
            NZ = LENR(I) - LENRL(I)
            IF (NZ.NE.0) GO TO 850
            LASTR(I) = 0
            GO TO 890
  850       IFIR = IFIRST(NZ)
            IFIRST(NZ) = I
            IF (IFIR) 860, 880, 870
  860       LASTR(I) = IFIR
            NEXTR(I) = 0
            GO TO 890
  870       LASTR(I) = LASTR(IFIR)
            NEXTR(I) = IFIR
            LASTR(IFIR) = I
            GO TO 890
  880       LASTR(I) = 0
            NEXTR(I) = 0
            NZMIN = MIN0(NZMIN,NZ)
  890     CONTINUE
c restore iq and nullify u part of old pivot row.
c    record the column permutation in lastc(jpiv) and the row
c    permutation in lastr(ipiv).
  900     IPC(JPIV) = -ISING
          LASTR(IPIV) = PIVOT
          IF (LENPIV.EQ.0) GO TO 980
          NZROW = NZROW - LENPIV
          JVAL = IJP1
          JZER = IPTR(IPIV)
          IPTR(IPIV) = 0
          DO 910 JCOUNT=1,LENPIV
            J = ICN(JVAL)
            IQ(J) = ICN(JZER)
            ICN(JZER) = 0
            JVAL = JVAL + 1
            JZER = JZER + 1
  910     CONTINUE
c adjust column ordering arrays.
          IF (NSRCH.GT.NN) GO TO 920
          DO 916 JJ=IJP1,PIVEND
            J = ICN(JJ)
            NZ = LENC(J)
            IF (NZ.NE.0) GO TO 914
            IPC(J) = 0
            GO TO 916
  914       NZMIN = MIN0(NZMIN,NZ)
  916     CONTINUE
          GO TO 980
  920     JJ = COLUPD
          DO 970 JDUMMY=1,NN
            J = JJ
            IF (J.EQ.NN+1) GO TO 980
            JJ = -NEXTC(J)
            NZ = LENC(J)
            IF (NZ.NE.0) GO TO 924
            IPC(J) = 0
            GO TO 970
  924       IFIR = IFIRST(NZ)
            LASTC(J) = 0
            IF (IFIR) 930, 940, 950
  930       IFIRST(NZ) = -J
            IFIR = -IFIR
            LASTC(IFIR) = J
            NEXTC(J) = IFIR
            GO TO 970
  940       IFIRST(NZ) = -J
            NEXTC(J) = 0
            GO TO 960
  950       LC = -LASTR(IFIR)
            LASTR(IFIR) = -J
            NEXTC(J) = LC
            IF (LC.NE.0) LASTC(LC) = J
  960       NZMIN = MIN0(NZMIN,NZ)
  970     CONTINUE
  980   CONTINUE
c ********************************************
c ****    end of main elimination loop    ****
c ********************************************
c
c reset iactiv to point to the beginning of the next block.
  990   IF (ILAST.NE.NN) IACTIV = IPTR(ILAST+1)
 1000 CONTINUE
c
c ********************************************
c ****    end of deomposition of block    ****
c ********************************************
c
c record singularity (if any) in iq array.
      IF (IRANK.EQ.NN) GO TO 1020
      DO 1010 I=1,NN
        IF (IPC(I).LT.0) GO TO 1010
        ISING = IPC(I)
        IQ(ISING) = -IQ(ISING)
        IPC(I) = -ISING
 1010 CONTINUE
c
c run through lu decomposition changing column indices to that of new
c     order and permuting lenr and lenrl arrays according to pivot
c     permutations.
 1020 ISTART = IDISP(1)
      IEND = IBEG - 1
      IF (IEND.LT.ISTART) GO TO 1040
      DO 1030 JJ=ISTART,IEND
        JOLD = ICN(JJ)
        ICN(JJ) = -IPC(JOLD)
 1030 CONTINUE
 1040 DO 1050 II=1,NN
        I = LASTR(II)
        NEXTR(I) = LENR(II)
        IPTR(I) = LENRL(II)
 1050 CONTINUE
      DO 1060 I=1,NN
        LENRL(I) = IPTR(I)
        LENR(I) = NEXTR(I)
 1060 CONTINUE
c
c update permutation arrays ip and iq.
      DO 1070 II=1,NN
        I = LASTR(II)
        J = -IPC(II)
        NEXTR(I) = IABS(IP(II)+0)
        IPTR(J) = IABS(IQ(II)+0)
 1070 CONTINUE
      DO 1080 I=1,NN
        IF (IP(I).LT.0) NEXTR(I) = -NEXTR(I)
        IP(I) = NEXTR(I)
        IF (IQ(I).LT.0) IPTR(I) = -IPTR(I)
        IQ(I) = IPTR(I)
 1080 CONTINUE
      IP(NN) = IABS(IP(NN)+0)
      IDISP(2) = IEND
      GO TO 1120
c
c   ***    error returns    ***
 1090 IDISP(2) = IACTIV
      IF (LP.EQ.0) GO TO 1120
      WRITE (LP,99996)
      GO TO 1110
 1100 IF (IFLAG.EQ.-5) IFLAG = -6
      IF (IFLAG.NE.-6) IFLAG = -3
      IDISP(2) = IACTIV
      IF (LP.EQ.0) GO TO 1120
      IF (IFLAG.EQ.-3) WRITE (LP,99995)
      IF (IFLAG.EQ.-6) WRITE (LP,99994)
 1110 PIVOT = PIVOT - ISTART + 1
      WRITE (LP,99993) PIVOT, NBLOCK, ISTART, ILAST
      IF (PIVOT.EQ.0) WRITE (LP,99992) MINIRN
c
c
 1120 RETURN
99999 FORMAT (54H ERROR RETURN FROM MA30A/AD BECAUSE MATRIX IS STRUCTUR,
     * 13HALLY SINGULAR)
99998 FORMAT (54H ERROR RETURN FROM MA30A/AD BECAUSE MATRIX IS NUMERICA,
     * 12HLLY SINGULAR)
99997 FORMAT (48H LU DECOMPOSITION DESTROYED TO CREATE MORE SPACE)
99996 FORMAT (54H ERROR RETURN FROM MA30A/AD BECAUSE LICN NOT BIG ENOUG,
     * 1HH)
99995 FORMAT (54H ERROR RETURN FROM MA30A/AD BECAUSE LIRN NOT BIG ENOUG,
     * 1HH)
99994 FORMAT (51H ERROR RETURN FROM MA30A/AD LIRN AND LICN TOO SMALL)
99993 FORMAT (10H AT STAGE , I5, 10H IN BLOCK , I5, 16H WITH FIRST ROW ,
     * I5, 14H AND LAST ROW , I5)
99992 FORMAT (34H TO CONTINUE SET LIRN TO AT LEAST , I8)
      END
c
c=====================================================================
c
      SUBROUTINE MA30BD(N, ICN, A, LICN, LENR, LENRL, IDISP, IP, IQ, W,
     * IW, IFLAG)
c ma30b/bd performs the lu decomposition of the diagonal blocks of a
c     new matrix paq of the same sparsity pattern, using information
c     from a previous call to ma30a/ad. the entries of the input
c     matrix  must already be in their final positions in the lu
c     decomposition structure.  this routine executes about five times
c     faster than ma30a/ad.
c
c we now describe the argument list for ma30b/bd. consult ma30a/ad for
c     further information on these parameters.
c n  is an integer variable set to the order of the matrix.
c icn is an integer array of length licn. it should be unchanged
c     since the last call to ma30a/ad. it is not altered by ma30b/bd.
c a  is a real/double precision array of length licn the user must set
c     entries idisp(1) to idisp(2) to contain the entries in the
c     diagonal blocks of the matrix paq whose column numbers are held
c     in icn, using corresponding positions. note that some zeros may
c     need to be held explicitly. on output entries idisp(1) to
c     idisp(2) of array a contain the lu decomposition of the diagonal
c     blocks of paq. entries a(1) to a(idisp(1)-1) are neither
c     required nor altered by ma30b/bd.
c licn  is an integer variable which must be set by the user to the
c     length of arrays a and icn. it is not altered by ma30b/bd.
c lenr,lenrl are integer arrays of length n. they should be
c     unchanged since the last call to ma30a/ad. they are not altered
c     by ma30b/bd.
c idisp  is an integer array of length 2. it should be unchanged since
c     the last call to ma30a/ad. it is not altered by ma30b/bd.
c ip,iq  are integer arrays of length n. they should be unchanged
c     since the last call to ma30a/ad. they are not altered by
c     ma30b/bd.
c w  is a real/double precision array of length n which is used as
c     workspace by ma30b/bd.
c iw  is an integer array of length n which is used as workspace by
c     ma30b/bd.
c iflag  is an integer variable. on output from ma30b/bd, iflag has
c     the value zero if the factorization was successful, has the
c     value i if pivot i was very small and has the value -i if an
c     unexpected singularity was detected at stage i of the
c     decomposition.
c
      DOUBLE PRECISION A(LICN), W(N), AU, EPS, ROWMAX, ZERO, ONE, RMIN,
     * TOL, BIG
      LOGICAL ABORT1, ABORT2, ABORT3, STAB, LBIG
      INTEGER IW(N), IDISP(2), PIVPOS
      INTEGER ICN(LICN), LENR(N), LENRL(N), IP(N), IQ(N)
c see block data for comments on variables in common.
      COMMON /MA30ED/ LP, ABORT1, ABORT2, ABORT3
      COMMON /MA30ID/ TOL, BIG, NDROP, NSRCH, LBIG
      COMMON /MA30GD/ EPS, RMIN
      DATA ZERO /0.0D0/, ONE /1.0D0/
      STAB = EPS.LE.ONE
      RMIN = EPS
      ISING = 0
      IFLAG = 0
      DO 10 I=1,N
        W(I) = ZERO
   10 CONTINUE
c set up pointers to the beginning of the rows.
      IW(1) = IDISP(1)
      IF (N.EQ.1) GO TO 25
      DO 20 I=2,N
        IW(I) = IW(I-1) + LENR(I-1)
   20 CONTINUE
c
c   ****   start  of main loop    ****
c at step i, row i of a is transformed to row i of l/u by adding
c     appropriate multiples of rows 1 to i-1.
c     .... using row-gauss elimination.
   25 DO 160 I=1,N
c istart is beginning of row i of a and row i of l.
        ISTART = IW(I)
c ifin is end of row i of a and row i of u.
        IFIN = ISTART + LENR(I) - 1
c ilend is end of row i of l.
        ILEND = ISTART + LENRL(I) - 1
        IF (ISTART.GT.ILEND) GO TO 90
c load row i of a into vector w.
        DO 30 JJ=ISTART,IFIN
          J = ICN(JJ)
          W(J) = A(JJ)
   30   CONTINUE
c
c add multiples of appropriate rows of  i to i-1  to row i.
        DO 70 JJ=ISTART,ILEND
          J = ICN(JJ)
c ipivj is position of pivot in row j.
          IPIVJ = IW(J) + LENRL(J)
c form multiplier au.
          AU = -W(J)/A(IPIVJ)
          IF (LBIG) BIG = DMAX1(DABS(AU),BIG)
          W(J) = AU
c au * row j (u part) is added to row i.
          IPIVJ = IPIVJ + 1
          JFIN = IW(J) + LENR(J) - 1
          IF (IPIVJ.GT.JFIN) GO TO 70
c innermost loop.
          IF (LBIG) GO TO 50
          DO 40 JAYJAY=IPIVJ,JFIN
            JAY = ICN(JAYJAY)
            W(JAY) = W(JAY) + AU*A(JAYJAY)
   40     CONTINUE
          GO TO 70
   50     DO 60 JAYJAY=IPIVJ,JFIN
            JAY = ICN(JAYJAY)
            W(JAY) = W(JAY) + AU*A(JAYJAY)
            BIG = DMAX1(DABS(W(JAY)),BIG)
   60     CONTINUE
   70   CONTINUE
c
c reload w back into a (now l/u)
        DO 80 JJ=ISTART,IFIN
          J = ICN(JJ)
          A(JJ) = W(J)
          W(J) = ZERO
   80   CONTINUE
c we now perform the stability checks.
   90   PIVPOS = ILEND + 1
        IF (IQ(I).GT.0) GO TO 140
c matrix had singularity at this point in ma30a/ad.
c is it the first such pivot in current block ?
        IF (ISING.EQ.0) ISING = I
c does current matrix have a singularity in the same place ?
        IF (PIVPOS.GT.IFIN) GO TO 100
        IF (A(PIVPOS).NE.ZERO) GO TO 170
c it does .. so set ising if it is not the end of the current block
c check to see that appropriate part of l/u is zero or null.
  100   IF (ISTART.GT.IFIN) GO TO 120
        DO 110 JJ=ISTART,IFIN
          IF (ICN(JJ).LT.ISING) GO TO 110
          IF (A(JJ).NE.ZERO) GO TO 170
  110   CONTINUE
  120   IF (PIVPOS.LE.IFIN) A(PIVPOS) = ONE
        IF (IP(I).GT.0 .AND. I.NE.N) GO TO 160
c end of current block ... reset zero pivots and ising.
        DO 130 J=ISING,I
          IF ((LENR(J)-LENRL(J)).EQ.0) GO TO 130
          JJ = IW(J) + LENRL(J)
          A(JJ) = ZERO
  130   CONTINUE
        ISING = 0
        GO TO 160
c matrix had non-zero pivot in ma30a/ad at this stage.
  140   IF (PIVPOS.GT.IFIN) GO TO 170
        IF (A(PIVPOS).EQ.ZERO) GO TO 170
        IF (.NOT.STAB) GO TO 160
        ROWMAX = ZERO
        DO 150 JJ=PIVPOS,IFIN
          ROWMAX = DMAX1(ROWMAX,DABS(A(JJ)))
  150   CONTINUE
        IF (DABS(A(PIVPOS))/ROWMAX.GE.RMIN) GO TO 160
        IFLAG = I
        RMIN = DABS(A(PIVPOS))/ROWMAX
c   ****    end of main loop    ****
  160 CONTINUE
c
      GO TO 180
c   ***   error return   ***
  170 IF (LP.NE.0) WRITE (LP,99999) I
      IFLAG = -I
c
  180 RETURN
99999 FORMAT (54H ERROR RETURN FROM MA30B/BD SINGULARITY DETECTED IN RO,
     * 1HW, I8)
      END
c
c======================================================================
c
      SUBROUTINE MA30CD(N, ICN, A, LICN, LENR, LENRL, LENOFF, IDISP, IP,
     * IQ, X, W, MTYPE)
c ma30c/cd uses the factors produced by ma30a/ad or ma30b/bd to solve
c     ax=b or a transpose x=b when the matrix p1*a*q1 (paq) is block
c     lower triangular (including the case of only one diagonal
c     block).
c
c we now describe the argument list for ma30c/cd.
c n  is an integer variable set to the order of the matrix. it is not
c     altered by the subroutine.
c icn is an integer array of length licn. entries idisp(1) to
c     idisp(2) should be unchanged since the last call to ma30a/ad. if
c     the matrix has more than one diagonal block, then column indices
c     corresponding to non-zeros in sub-diagonal blocks of paq must
c     appear in positions 1 to idisp(1)-1. for the same row those
c     entries must be contiguous, with those in row i preceding those
c     in row i+1 (i=1,...,n-1) and no wasted space between rows.
c     entries may be in any order within each row. it is not altered
c     by ma30c/cd.
c a  is a real/double precision array of length licn.  entries
c     idisp(1) to idisp(2) should be unchanged since the last call to
c     ma30a/ad or ma30b/bd.  if the matrix has more than one diagonal
c     block, then the values of the non-zeros in sub-diagonal blocks
c     must be in positions 1 to idisp(1)-1 in the order given by icn.
c     it is not altered by ma30c/cd.
c licn  is an integer variable set to the size of arrays icn and a.
c     it is not altered by ma30c/cd.
c lenr,lenrl are integer arrays of length n which should be
c     unchanged since the last call to ma30a/ad. they are not altered
c     by ma30c/cd.
c lenoff  is an integer array of length n. if the matrix paq (or
c     p1*a*q1) has more than one diagonal block, then lenoff(i),
c     i=1,...,n should be set to the number of non-zeros in row i of
c     the matrix paq which are in sub-diagonal blocks.  if there is
c     only one diagonal block then lenoff(1) may be set to -1, in
c     which case the other entries of lenoff are never accessed. it is
c     not altered by ma30c/cd.
c idisp  is an integer array of length 2 which should be unchanged
c     since the last call to ma30a/ad. it is not altered by ma30c/cd.
c ip,iq are integer arrays of length n which should be unchanged
c     since the last call to ma30a/ad. they are not altered by
c     ma30c/cd.
c x is a real/double precision array of length n. it must be set by
c     the user to the values of the right hand side vector b for the
c     equations being solved.  on exit from ma30c/cd it will be equal
c     to the solution x required.
c w  is a real/double precision array of length n which is used as
c     workspace by ma30c/cd.
c mtype is an integer variable which must be set by the user. if
c     mtype=1, then the solution to the system ax=b is returned; any
c     other value for mtype will return the solution to the system a
c     transpose x=b. it is not altered by ma30c/cd.
c
      DOUBLE PRECISION A(LICN), X(N), W(N), WII, WI, RESID, ZERO
      LOGICAL NEG, NOBLOC
      INTEGER IDISP(2)
      INTEGER ICN(LICN), LENR(N), LENRL(N), LENOFF(N), IP(N), IQ(N)
c see block data for comments on variables in common.
      COMMON /MA30HD/ RESID
      DATA ZERO /0.0D0/
c
c the final value of resid is the maximum residual for an inconsistent
c     set of equations.
      RESID = ZERO
c nobloc is .true. if subroutine block has been used previously and
c     is .false. otherwise.  the value .false. means that lenoff
c     will not be subsequently accessed.
      NOBLOC = LENOFF(1).LT.0
      IF (MTYPE.NE.1) GO TO 140
c
c we now solve   a * x = b.
c neg is used to indicate when the last row in a block has been
c     reached.  it is then set to true whereafter backsubstitution is
c     performed on the block.
      NEG = .FALSE.
c ip(n) is negated so that the last row of the last block can be
c     recognised.  it is reset to its positive value on exit.
      IP(N) = -IP(N)
c preorder vector ... w(i) = x(ip(i))
      DO 10 II=1,N
        I = IP(II)
        I = IABS(I)
        W(II) = X(I)
   10 CONTINUE
c lt holds the position of the first non-zero in the current row of the
c     off-diagonal blocks.
      LT = 1
c ifirst holds the index of the first row in the current block.
      IFIRST = 1
c iblock holds the position of the first non-zero in the current row
c     of the lu decomposition of the diagonal blocks.
      IBLOCK = IDISP(1)
c if i is not the last row of a block, then a pass through this loop
c     adds the inner product of row i of the off-diagonal blocks and w
c     to w and performs forward elimination using row i of the lu
c     decomposition.   if i is the last row of a block then, after
c     performing these aforementioned operations, backsubstitution is
c     performed using the rows of the block.
      DO 120 I=1,N
        WI = W(I)
        IF (NOBLOC) GO TO 30
        IF (LENOFF(I).EQ.0) GO TO 30
c operations using lower triangular blocks.
c ltend is the end of row i in the off-diagonal blocks.
        LTEND = LT + LENOFF(I) - 1
        DO 20 JJ=LT,LTEND
          J = ICN(JJ)
          WI = WI - A(JJ)*W(J)
   20   CONTINUE
c lt is set the beginning of the next off-diagonal row.
        LT = LTEND + 1
c set neg to .true. if we are on the last row of the block.
   30   IF (IP(I).LT.0) NEG = .TRUE.
        IF (LENRL(I).EQ.0) GO TO 50
c forward elimination phase.
c iend is the end of the l part of row i in the lu decomposition.
        IEND = IBLOCK + LENRL(I) - 1
        DO 40 JJ=IBLOCK,IEND
          J = ICN(JJ)
          WI = WI + A(JJ)*W(J)
   40   CONTINUE
c iblock is adjusted to point to the start of the next row.
   50   IBLOCK = IBLOCK + LENR(I)
        W(I) = WI
        IF (.NOT.NEG) GO TO 120
c back substitution phase.
c j1 is position in a/icn after end of block beginning in row ifirst
c     and ending in row i.
        J1 = IBLOCK
c are there any singularities in this block?  if not, continue with
c     the backsubstitution.
        IB = I
        IF (IQ(I).GT.0) GO TO 70
        DO 60 III=IFIRST,I
          IB = I - III + IFIRST
          IF (IQ(IB).GT.0) GO TO 70
          J1 = J1 - LENR(IB)
          RESID = DMAX1(RESID,DABS(W(IB)))
          W(IB) = ZERO
   60   CONTINUE
c entire block is singular.
        GO TO 110
c each pass through this loop performs the back-substitution
c     operations for a single row, starting at the end of the block and
c     working through it in reverse order.
   70   DO 100 III=IFIRST,IB
          II = IB - III + IFIRST
c j2 is end of row ii.
          J2 = J1 - 1
c j1 is beginning of row ii.
          J1 = J1 - LENR(II)
c jpiv is the position of the pivot in row ii.
          JPIV = J1 + LENRL(II)
          JPIVP1 = JPIV + 1
c jump if row  ii of u has no non-zeros.
          IF (J2.LT.JPIVP1) GO TO 90
          WII = W(II)
          DO 80 JJ=JPIVP1,J2
            J = ICN(JJ)
            WII = WII - A(JJ)*W(J)
   80     CONTINUE
          W(II) = WII
   90     W(II) = W(II)/A(JPIV)
  100   CONTINUE
  110   IFIRST = I + 1
        NEG = .FALSE.
  120 CONTINUE
c
c reorder solution vector ... x(i) = w(iqinverse(i))
      DO 130 II=1,N
        I = IQ(II)
        I = IABS(I)
        X(I) = W(II)
  130 CONTINUE
      IP(N) = -IP(N)
      GO TO 320
c
c
c we now solve   atranspose * x = b.
c preorder vector ... w(i)=x(iq(i))
  140 DO 150 II=1,N
        I = IQ(II)
        I = IABS(I)
        W(II) = X(I)
  150 CONTINUE
c lj1 points to the beginning the current row in the off-diagonal
c     blocks.
      LJ1 = IDISP(1)
c iblock is initialized to point to the beginning of the block after
c     the last one ]
      IBLOCK = IDISP(2) + 1
c ilast is the last row in the current block.
      ILAST = N
c iblend points to the position after the last non-zero in the
c     current block.
      IBLEND = IBLOCK
c each pass through this loop operates with one diagonal block and
c     the off-diagonal part of the matrix corresponding to the rows
c     of this block.  the blocks are taken in reverse order and the
c     number of times the loop is entered is min(n,no. blocks+1).
      DO 290 NUMBLK=1,N
        IF (ILAST.EQ.0) GO TO 300
        IBLOCK = IBLOCK - LENR(ILAST)
c this loop finds the index of the first row in the current block..
c     it is first and iblock is set to the position of the beginning
c     of this first row.
        DO 160 K=1,N
          II = ILAST - K
          IF (II.EQ.0) GO TO 170
          IF (IP(II).LT.0) GO TO 170
          IBLOCK = IBLOCK - LENR(II)
  160   CONTINUE
  170   IFIRST = II + 1
c j1 points to the position of the beginning of row i (lt part) or pivot
        J1 = IBLOCK
c forward elimination.
c each pass through this loop performs the operations for one row of the
c     block.  if the corresponding entry of w is zero then the
c     operations can be avoided.
        DO 210 I=IFIRST,ILAST
          IF (W(I).EQ.ZERO) GO TO 200
c jump if row i singular.
          IF (IQ(I).LT.0) GO TO 220
c j2 first points to the pivot in row i and then is made to point to the
c     first non-zero in the u transpose part of the row.
          J2 = J1 + LENRL(I)
          WI = W(I)/A(J2)
          IF (LENR(I)-LENRL(I).EQ.1) GO TO 190
          J2 = J2 + 1
c j3 points to the end of row i.
          J3 = J1 + LENR(I) - 1
          DO 180 JJ=J2,J3
            J = ICN(JJ)
            W(J) = W(J) - A(JJ)*WI
  180     CONTINUE
  190     W(I) = WI
  200     J1 = J1 + LENR(I)
  210   CONTINUE
        GO TO 240
c deals with rest of block which is singular.
  220   DO 230 II=I,ILAST
          RESID = DMAX1(RESID,DABS(W(II)))
          W(II) = ZERO
  230   CONTINUE
c back substitution.
c this loop does the back substitution on the rows of the block in
c     the reverse order doing it simultaneously on the l transpose part
c     of the diagonal blocks and the off-diagonal blocks.
  240   J1 = IBLEND
        DO 280 IBACK=IFIRST,ILAST
          I = ILAST - IBACK + IFIRST
c j1 points to the beginning of row i.
          J1 = J1 - LENR(I)
          IF (LENRL(I).EQ.0) GO TO 260
c j2 points to the end of the l transpose part of row i.
          J2 = J1 + LENRL(I) - 1
          DO 250 JJ=J1,J2
            J = ICN(JJ)
            W(J) = W(J) + A(JJ)*W(I)
  250     CONTINUE
  260     IF (NOBLOC) GO TO 280
c operations using lower triangular blocks.
          IF (LENOFF(I).EQ.0) GO TO 280
c lj2 points to the end of row i of the off-diagonal blocks.
          LJ2 = LJ1 - 1
c lj1 points to the beginning of row i of the off-diagonal blocks.
          LJ1 = LJ1 - LENOFF(I)
          DO 270 JJ=LJ1,LJ2
            J = ICN(JJ)
            W(J) = W(J) - A(JJ)*W(I)
  270     CONTINUE
  280   CONTINUE
        IBLEND = J1
        ILAST = IFIRST - 1
  290 CONTINUE
c reorder solution vector ... x(i)=w(ipinverse(i))
  300 DO 310 II=1,N
        I = IP(II)
        I = IABS(I)
        X(I) = W(II)
  310 CONTINUE
c
  320 RETURN
      END
      SUBROUTINE MA30DD(A, ICN, IPTR, N, IACTIV, ITOP, REALS)
c this subroutine performs garbage collection operations on the
c     arrays a, icn and irn.
c iactiv is the first position in arrays a/icn from which the compress
c     starts.  on exit, iactiv equals the position of the first entry
c     in the compressed part of a/icn
c
      DOUBLE PRECISION A(ITOP)
      LOGICAL REALS
      INTEGER IPTR(N)
      INTEGER ICN(ITOP)
c see block data for comments on variables in common.
      COMMON /MA30FD/ IRNCP, ICNCP, IRANK, MINIRN, MINICN
c
      IF (REALS) ICNCP = ICNCP + 1
      IF (.NOT.REALS) IRNCP = IRNCP + 1
c set the first non-zero entry in each row to the negative of the
c     row/col number and hold this row/col index in the row/col
c     pointer.  this is so that the beginning of each row/col can
c     be recognized in the subsequent scan.
      DO 10 J=1,N
        K = IPTR(J)
        IF (K.LT.IACTIV) GO TO 10
        IPTR(J) = ICN(K)
        ICN(K) = -J
   10 CONTINUE
      KN = ITOP + 1
      KL = ITOP - IACTIV + 1
c go through arrays in reverse order compressing to the back so
c     that there are no zeros held in positions iactiv to itop in icn.
c     reset first entry of each row/col and pointer array iptr.
      DO 30 K=1,KL
        JPOS = ITOP - K + 1
        IF (ICN(JPOS).EQ.0) GO TO 30
        KN = KN - 1
        IF (REALS) A(KN) = A(JPOS)
        IF (ICN(JPOS).GE.0) GO TO 20
c first non-zero of row/col has been located
        J = -ICN(JPOS)
        ICN(JPOS) = IPTR(J)
        IPTR(J) = KN
   20   ICN(KN) = ICN(JPOS)
   30 CONTINUE
      IACTIV = KN
      RETURN
      END
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc13d
      SUBROUTINE MC13D(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW)
      INTEGER IP(N)
      INTEGER ICN(LICN),LENR(N),IOR(N),IB(N),IW(N,3)
      CALL MC13E(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW(1,1),IW(1,2),IW(1,3))
      RETURN
      END
      SUBROUTINE MC13E(N,ICN,LICN,IP,LENR,ARP,IB,NUM,LOWL,NUMB,PREV)
      INTEGER STP,DUMMY
      INTEGER IP(N)
c
c arp(i) is one less than the number of unsearched edges leaving
c     node i.  at the end of the algorithm it is set to a
c     permutation which puts the matrix in block lower
c     triangular form.
c ib(i) is the position in the ordering of the start of the ith
c     block.  ib(n+1-i) holds the node number of the ith node
c     on the stack.
c lowl(i) is the smallest stack position of any node to which a path
c     from node i has been found.  it is set to n+1 when node i
c     is removed from the stack.
c numb(i) is the position of node i in the stack if it is on
c     it, is the permuted order of node i for those nodes
c     whose final position has been found and is otherwise zero.
c prev(i) is the node at the end of the path when node i was
c     placed on the stack.
      INTEGER ICN(LICN),LENR(N),ARP(N),IB(N),LOWL(N),NUMB(N),
     1PREV(N)
c
c
c   icnt is the number of nodes whose positions in final ordering have
c     been found.
      ICNT=0
c num is the number of blocks that have been found.
      NUM=0
      NNM1=N+N-1
c
c initialization of arrays.
      DO 20 J=1,N
      NUMB(J)=0
      ARP(J)=LENR(J)-1
   20 CONTINUE
c
c
      DO 120 ISN=1,N
c look for a starting node
      IF (NUMB(ISN).NE.0) GO TO 120
      IV=ISN
c ist is the number of nodes on the stack ... it is the stack pointer.
      IST=1
c put node iv at beginning of stack.
      LOWL(IV)=1
      NUMB(IV)=1
      IB(N)=IV
c
c the body of this loop puts a new node on the stack or backtracks.
      DO 110 DUMMY=1,NNM1
      I1=ARP(IV)
c have all edges leaving node iv been searched.
      IF (I1.LT.0) GO TO 60
      I2=IP(IV)+LENR(IV)-1
      I1=I2-I1
c
c look at edges leaving node iv until one enters a new node or
c     all edges are exhausted.
      DO 50 II=I1,I2
      IW=ICN(II)
c has node iw been on stack already.
      IF (NUMB(IW).EQ.0) GO TO 100
c update value of lowl(iv) if necessary.
  50  LOWL(IV)=MIN0(LOWL(IV),LOWL(IW))
c
c there are no more edges leaving node iv.
      ARP(IV)=-1
c is node iv the root of a block.
   60 IF (LOWL(IV).LT.NUMB(IV)) GO TO 90
c
c order nodes in a block.
      NUM=NUM+1
      IST1=N+1-IST
      LCNT=ICNT+1
c peel block off the top of the stack starting at the top and
c     working down to the root of the block.
      DO 70 STP=IST1,N
      IW=IB(STP)
      LOWL(IW)=N+1
      ICNT=ICNT+1
      NUMB(IW)=ICNT
      IF (IW.EQ.IV) GO TO 80
   70 CONTINUE
   80 IST=N-STP
      IB(NUM)=LCNT
c are there any nodes left on the stack.
      IF (IST.NE.0) GO TO 90
c have all the nodes been ordered.
      IF (ICNT.LT.N) GO TO 120
      GO TO 130
c
c backtrack to previous node on path.
   90 IW=IV
      IV=PREV(IV)
c update value of lowl(iv) if necessary.
      LOWL(IV)=MIN0(LOWL(IV),LOWL(IW))
      GO TO 110
c
c put new node on the stack.
 100  ARP(IV)=I2-II-1
      PREV(IW)=IV
      IV=IW
      IST=IST+1
      LOWL(IV)=IST
      NUMB(IV)=IST
      K=N+1-IST
      IB(K)=IV
  110 CONTINUE
c
  120 CONTINUE
c
c
c put permutation in the required form.
  130 DO 140 I=1,N
      II=NUMB(I)
 140  ARP(II)=I
      RETURN
      END
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc20ad mc20bd
      SUBROUTINE MC20AD(NC,MAXA,A,INUM,JPTR,JNUM,JDISP)
c
      INTEGER   INUM(MAXA),JNUM(MAXA)
      DOUBLE PRECISION A(MAXA),ACE,ACEP
      DIMENSION JPTR(NC)
c
c     ******************************************************************
c
      NULL=-JDISP
c**      clear jptr
      DO 60 J=1,NC
   60 JPTR(J)=0
c**      count the number of elements in each column.
      DO 120 K=1,MAXA
      J=JNUM(K)+JDISP
      JPTR(J)=JPTR(J)+1
  120 CONTINUE
c**      set the jptr array
      K=1
      DO 150 J=1,NC
      KR=K+JPTR(J)
      JPTR(J)=K
  150 K=KR
c
c**      reorder the elements into column order.  the algorithm is an
c        in-place sort and is of order maxa.
      DO 230 I=1,MAXA
c        establish the current entry.
      JCE=JNUM(I)+JDISP
      IF(JCE.EQ.0) GO TO 230
      ACE=A(I)
      ICE=INUM(I)
c        clear the location vacated.
      JNUM(I)=NULL
c        chain from current entry to store items.
      DO 200 J=1,MAXA
c        current entry not in correct position.  determine correct
c        position to store entry.
      LOC=JPTR(JCE)
      JPTR(JCE)=JPTR(JCE)+1
c        save contents of that location.
      ACEP=A(LOC)
      ICEP=INUM(LOC)
      JCEP=JNUM(LOC)
c        store current entry.
      A(LOC)=ACE
      INUM(LOC)=ICE
      JNUM(LOC)=NULL
c        check if next current entry needs to be processed.
      IF(JCEP.EQ.NULL) GO TO 230
c        it does.  copy into current entry.
      ACE=ACEP
      ICE=ICEP
  200 JCE=JCEP+JDISP
c
  230 CONTINUE
c
c**      reset jptr vector.
      JA=1
      DO 250 J=1,NC
      JB=JPTR(J)
      JPTR(J)=JA
  250 JA=JB
      RETURN
      END
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc21a
      SUBROUTINE MC21A(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
      INTEGER IP(N)
      INTEGER ICN(LICN),LENR(N),IPERM(N),IW(N,4)
      CALL MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),IW(1,3),
     1IW(1,4))
      RETURN
      END
      SUBROUTINE MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
      INTEGER IP(N)
c   pr(i) is the previous row to i in the depth first search.
c it is used as a work array in the sorting algorithm.
c   elements (iperm(i),i) i=1, ... n  are non-zero at the end of the
c algorithm unless n assignments have not been made.  in which case
c (iperm(i),i) will be zero for n-numnz entries.
c   cv(i) is the most recent row extension at which column i
c was visited.
c   arp(i) is one less than the number of non-zeros in row i
c which have not been scanned when looking for a cheap assignment.
c   out(i) is one less than the number of non-zeros in row i
c which have not been scanned during one pass through the main loop.
      INTEGER ICN(LICN),LENR(N),IPERM(N),PR(N),CV(N),
     1ARP(N),OUT(N)
c
c   initialization of arrays.
      DO 10 I=1,N
      ARP(I)=LENR(I)-1
      CV(I)=0
   10 IPERM(I)=0
      NUMNZ=0
c
c
c   main loop.
c   each pass round this loop either results in a new assignment
c or gives a row with no assignment.
      DO 130 JORD=1,N
      J=JORD
      PR(J)=-1
      DO 100 K=1,JORD
c look for a cheap assignment
      IN1=ARP(J)
      IF (IN1.LT.0) GO TO 60
      IN2=IP(J)+LENR(J)-1
      IN1=IN2-IN1
      DO 50 II=IN1,IN2
      I=ICN(II)
      IF (IPERM(I).EQ.0) GO TO 110
   50 CONTINUE
c   no cheap assignment in row.
      ARP(J)=-1
c   begin looking for assignment chain starting with row j.
   60 OUT(J)=LENR(J)-1
c inner loop.  extends chain by one or backtracks.
      DO 90 KK=1,JORD
      IN1=OUT(J)
      IF (IN1.LT.0) GO TO 80
      IN2=IP(J)+LENR(J)-1
      IN1=IN2-IN1
c forward scan.
      DO 70 II=IN1,IN2
      I=ICN(II)
      IF (CV(I).EQ.JORD) GO TO 70
c   column i has not yet been accessed during this pass.
      J1=J
      J=IPERM(I)
      CV(I)=JORD
      PR(J)=J1
      OUT(J1)=IN2-II-1
      GO TO 100
   70 CONTINUE
c
c   backtracking step.
   80 J=PR(J)
      IF (J.EQ.-1) GO TO 130
   90 CONTINUE
c
  100 CONTINUE
c
c   new assignment is made.
  110 IPERM(I)=J
      ARP(J)=IN2-II-1
      NUMNZ=NUMNZ+1
      DO 120 K=1,JORD
      J=PR(J)
      IF (J.EQ.-1) GO TO 130
      II=IP(J)+LENR(J)-OUT(J)-2
      I=ICN(II)
      IPERM(I)=J
  120 CONTINUE
c
  130 CONTINUE
c
c   if matrix is structurally singular, we now complete the
c permutation iperm.
      IF (NUMNZ.EQ.N) RETURN
      DO 140 I=1,N
  140 ARP(I)=0
      K=0
      DO 160 I=1,N
      IF (IPERM(I).NE.0) GO TO 150
      K=K+1
      OUT(K)=I
      GO TO 160
  150 J=IPERM(I)
      ARP(J)=I
  160 CONTINUE
      K=0
      DO 170 I=1,N
      IF (ARP(I).NE.0) GO TO 170
      K=K+1
      IOUTK=OUT(K)
      IPERM(IOUTK)=I
  170 CONTINUE
      RETURN
      END
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc22ad
      SUBROUTINE MC22AD(N,ICN,A,NZ,LENROW,IP,IQ,IW,IW1)
      DOUBLE PRECISION A(NZ),AVAL
      INTEGER IW(N,2)
      INTEGER   ICN(NZ),LENROW(N),IP(N),IQ(N),IW1(NZ)
      IF (NZ.LE.0) GO TO 1000
      IF (N.LE.0) GO TO 1000
c set start of row i in iw(i,1) and lenrow(i) in iw(i,2)
      IW(1,1)=1
      IW(1,2)=LENROW(1)
      DO 10 I=2,N
      IW(I,1)=IW(I-1,1)+LENROW(I-1)
 10   IW(I,2)=LENROW(I)
c permute lenrow according to ip.  set off-sets for new position
c     of row iold in iw(iold,1) and put old row indices in iw1 in
c     positions corresponding to the new position of this row in a/icn.
      JJ=1
      DO 20 I=1,N
      IOLD=IP(I)
      IOLD=IABS(IOLD)
      LENGTH=IW(IOLD,2)
      LENROW(I)=LENGTH
      IF (LENGTH.EQ.0) GO TO 20
      IW(IOLD,1)=IW(IOLD,1)-JJ
      J2=JJ+LENGTH-1
      DO 15 J=JJ,J2
 15   IW1(J)=IOLD
      JJ=J2+1
 20   CONTINUE
c set inverse permutation to iq in iw(.,2).
      DO 30 I=1,N
      IOLD=IQ(I)
      IOLD=IABS(IOLD)
 30   IW(IOLD,2)=I
c permute a and icn in place, changing to new column numbers.
c
c ***   main loop   ***
c each pass through this loop places a closed chain of column indices
c     in their new (and final) positions ... this is recorded by
c     setting the iw1 entry to zero so that any which are subsequently
c     encountered during this major scan can be bypassed.
      DO 200 I=1,NZ
      IOLD=IW1(I)
      IF (IOLD.EQ.0) GO TO 200
      IPOS=I
      JVAL=ICN(I)
c if row iold is in same positions after permutation go to 150.
      IF (IW(IOLD,1).EQ.0) GO TO 150
      AVAL=A(I)
c **  chain loop  **
c each pass through this loop places one (permuted) column index
c     in its final position  .. viz. ipos.
      DO 100 ICHAIN=1,NZ
c newpos is the original position in a/icn of the element to be placed
c in position ipos.  it is also the position of the next element in
c     the chain.
      NEWPOS=IPOS+IW(IOLD,1)
c is chain complete ?
      IF (NEWPOS.EQ.I) GO TO 130
      A(IPOS)=A(NEWPOS)
      JNUM=ICN(NEWPOS)
      ICN(IPOS)=IW(JNUM,2)
      IPOS=NEWPOS
      IOLD=IW1(IPOS)
      IW1(IPOS)=0
c **  end of chain loop  **
 100  CONTINUE
 130  A(IPOS)=AVAL
 150  ICN(IPOS)=IW(JVAL,2)
c ***   end of main loop   ***
 200  CONTINUE
c
 1000 RETURN
      END
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc23ad
c###### calls   mc13    mc21
      SUBROUTINE MC23AD(N,ICN,A,LICN,LENR,IDISP,IP,IQ,LENOFF,IW,IW1)
      DOUBLE PRECISION A(LICN)
      INTEGER IDISP(2),IW1(N,2)
      LOGICAL ABORT
      INTEGER   ICN(LICN),LENR(N),IP(N),IQ(N),LENOFF(N),IW(N,5)
      COMMON /MC23BD/ LP,NUMNZ,NUM,LARGE,ABORT
c input ... n,icn .. a,icn,lenr ....
c
c set up pointers iw(.,1) to the beginning of the rows and set lenoff
c     equal to lenr.
      IW1(1,1)=1
      LENOFF(1)=LENR(1)
      IF (N.EQ.1) GO TO 20
      DO 10 I=2,N
      LENOFF(I)=LENR(I)
   10 IW1(I,1)=IW1(I-1,1)+LENR(I-1)
c idisp(1) points to the first position in a/icn after the
c     off-diagonal blocks and untreated rows.
   20 IDISP(1)=IW1(N,1)+LENR(N)
c
c find row permutation ip to make diagonal zero-free.
      CALL MC21A(N,ICN,LICN,IW1,LENR,IP,NUMNZ,IW)
c
c possible error return for structurally singular matrices.
      IF (NUMNZ.NE.N.AND.ABORT) GO TO 170
c
c iw1(.,2) and lenr are permutations of iw1(.,1) and lenr/lenoff
c     suitable for entry
c     to mc13d since matrix with these row pointer and length arrays
c     has maximum number of non-zeros on the diagonal.
      DO 30 II=1,N
      I=IP(II)
      IW1(II,2)=IW1(I,1)
   30 LENR(II)=LENOFF(I)
c
c find symmetric permutation iq to block lower triangular form.
      CALL MC13D(N,ICN,LICN,IW1(1,2),LENR,IQ,IW(1,4),NUM,IW)
c
      IF (NUM.NE.1) GO TO 60
c
c action taken if matrix is irreducible.
c whole matrix is just moved to the end of the storage.
      DO 40 I=1,N
      LENR(I)=LENOFF(I)
      IP(I)=I
   40 IQ(I)=I
      LENOFF(1)=-1
c idisp(1) is the first position after the last element in the
c     off-diagonal blocks and untreated rows.
      NZ=IDISP(1)-1
      IDISP(1)=1
c idisp(2) is the position in a/icn of the first element in the
c     diagonal blocks.
      IDISP(2)=LICN-NZ+1
      LARGE=N
      IF (NZ.EQ.LICN) GO TO 230
      DO 50 K=1,NZ
      J=NZ-K+1
      JJ=LICN-K+1
      A(JJ)=A(J)
   50 ICN(JJ)=ICN(J)
c 230 = return
      GO TO 230
c
c data structure reordered.
c
c form composite row permutation ... ip(i) = ip(iq(i)).
   60 DO 70 II=1,N
      I=IQ(II)
   70 IW(II,1)=IP(I)
      DO 80 I=1,N
   80 IP(I)=IW(I,1)
c
c run through blocks in reverse order separating diagonal blocks
c     which are moved to the end of the storage.  elements in
c     off-diagonal blocks are left in place unless a compress is
c     necessary.
c
c ibeg indicates the lowest value of j for which icn(j) has been
c     set to zero when element in position j was moved to the
c     diagonal block part of storage.
      IBEG=LICN+1
c iend is the position of the first element of those treated rows
c     which are in diagonal blocks.
      IEND=LICN+1
c large is the dimension of the largest block encountered so far.
      LARGE=0
c
c num is the number of diagonal blocks.
      DO 150 K=1,NUM
      IBLOCK=NUM-K+1
c i1 is first row (in permuted form) of block iblock.
c i2 is last row (in permuted form) of block iblock.
      I1=IW(IBLOCK,4)
      I2=N
      IF (K.NE.1) I2=IW(IBLOCK+1,4)-1
      LARGE=MAX0(LARGE,I2-I1+1)
c go through the rows of block iblock in the reverse order.
      DO 140 II=I1,I2
      INEW=I2-II+I1
c we now deal with row inew in permuted form (row iold in original
c     matrix).
      IOLD=IP(INEW)
c if there is space to move up diagonal block portion of row go to 110
      IF (IEND-IDISP(1).GE.LENOFF(IOLD)) GO TO 110
c
c in-line compress.
c moves separated off-diagonal elements and untreated rows to
c     front of storage.
      JNPOS=IBEG
      ILEND=IDISP(1)-1
      IF (ILEND.LT.IBEG) GO TO 190
      DO 90 J=IBEG,ILEND
      IF (ICN(J).EQ.0) GO TO 90
      ICN(JNPOS)=ICN(J)
      A(JNPOS)=A(J)
      JNPOS=JNPOS+1
   90 CONTINUE
      IDISP(1)=JNPOS
      IF (IEND-JNPOS.LT.LENOFF(IOLD)) GO TO 190
      IBEG=LICN+1
c reset pointers to the beginning of the rows.
      DO 100 I=2,N
  100 IW1(I,1)=IW1(I-1,1)+LENOFF(I-1)
c
c row iold is now split into diag. and off-diag. parts.
  110 IROWB=IW1(IOLD,1)
      LENI=0
      IROWE=IROWB+LENOFF(IOLD)-1
c backward scan of whole of row iold (in original matrix).
      IF (IROWE.LT.IROWB) GO TO 130
      DO 120 JJ=IROWB,IROWE
      J=IROWE-JJ+IROWB
      JOLD=ICN(J)
c iw(.,2) holds the inverse permutation to iq.
c     ..... it was set to this in mc13d.
      JNEW=IW(JOLD,2)
c if (jnew.lt.i1) then ....
c element is in off-diagonal block and so is left in situ.
      IF (JNEW.LT.I1) GO TO 120
c element is in diagonal block and is moved to the end of the storage.
      IEND=IEND-1
      A(IEND)=A(J)
      ICN(IEND)=JNEW
      IBEG=MIN0(IBEG,J)
      ICN(J)=0
      LENI=LENI+1
  120 CONTINUE
c
      LENOFF(IOLD)=LENOFF(IOLD)-LENI
  130 LENR(INEW)=LENI
  140 CONTINUE
c
      IP(I2)=-IP(I2)
  150 CONTINUE
c resets ip(n) to positive value.
      IP(N)=-IP(N)
c idisp(2) is position of first element in diagonal blocks.
      IDISP(2)=IEND
c
c this compress is used to move all off-diagonal elements to the
c     front of the storage.
      IF (IBEG.GT.LICN) GO TO 230
      JNPOS=IBEG
      ILEND=IDISP(1)-1
      DO 160 J=IBEG,ILEND
      IF (ICN(J).EQ.0) GO TO 160
      ICN(JNPOS)=ICN(J)
      A(JNPOS)=A(J)
      JNPOS=JNPOS+1
  160 CONTINUE
c idisp(1) is first position after last element of off-diagonal blocks.
      IDISP(1)=JNPOS
      GO TO 230
c
c
c error return
  170 IF (LP.NE.0) WRITE(LP,180) NUMNZ
  180 FORMAT(33X,41H MATRIX IS STRUCTURALLY SINGULAR, RANK = ,I6)
      IDISP(1)=-1
      GO TO 210
  190 IF (LP.NE.0) WRITE(LP,200) N
  200 FORMAT(33X,33H LICN NOT BIG ENOUGH INCREASE BY ,I6)
      IDISP(1)=-2
  210 IF (LP.NE.0) WRITE(LP,220)
  220 FORMAT(33H+ERROR RETURN FROM MC23AD BECAUSE)
c
  230 RETURN
      END
c######date   01 jan 1984     copyright ukaea, harwell.
c######alias mc24ad
      SUBROUTINE MC24AD(N,ICN,A,LICN,LENR,LENRL,W)
      DOUBLE PRECISION A(LICN),W(N),AMAXL,WROWL,AMAXU,ZERO
      INTEGER   ICN(LICN),LENR(N),LENRL(N)
      DATA ZERO/0.0D0/
      AMAXL=ZERO
      DO 10 I=1,N
 10   W(I)=ZERO
      J0=1
      DO 100 I=1,N
      IF (LENR(I).EQ.0) GO TO 100
      J2=J0+LENR(I)-1
      IF (LENRL(I).EQ.0) GO TO 50
c calculation of 1-norm of l.
      J1=J0+LENRL(I)-1
      WROWL=ZERO
      DO 30 JJ=J0,J1
 30   WROWL=WROWL+DABS(A(JJ))
c amaxl is the maximum norm of columns of l so far found.
      AMAXL=DMAX1(AMAXL,WROWL)
      J0=J1+1
c calculation of norms of columns of u (max-norms).
 50   J0=J0+1
      IF (J0.GT.J2) GO TO 90
      DO 80 JJ=J0,J2
      J=ICN(JJ)
 80   W(J)=DMAX1(DABS(A(JJ)),W(J))
 90   J0=J2+1
 100  CONTINUE
c amaxu is set to maximum max-norm of columns of u.
      AMAXU=ZERO
      DO 200 I=1,N
 200  AMAXU=DMAX1(AMAXU,W(I))
c grofac is max u max-norm times max l 1-norm.
      W(1)=AMAXL*AMAXU
      RETURN
      END
      BLOCK DATA MABLD1
c
c comments on all the common block variables are given here even
c     though some are not initialized by block data.
c lp,mp are used by the subroutine as the unit numbers for its warning
c     and diagnostic messages. default value for both is 6 (for line
c     printer output). the user can either reset them to a different
c     stream number or suppress the output by setting them to zero.
c     while lp directs the output of error diagnostics from the
c     principal subroutines and internally called subroutines, mp
c     controls only the output of a message which warns the user that he
c     has input two or more non-zeros a(i), . . ,a(k) with the same row
c     and column indices.  the action taken in this case is to proceed
c     using a numerical value of a(i)+...+a(k). in the absence of other
c     errors, iflag will equal -14 on exit.
c lblock is a logical variable which controls an option of first
c     preordering the matrix to block lower triangular form (using
c     harwell subroutine mc23a). the preordering is performed if lblock
c     is equal to its default value of .true. if lblock is set to
c     .false. , the option is not invoked and the space allocated to
c     ikeep can be reduced to 4*n+1.
c grow is a logical variable. if it is left at its default value of
c     .true. , then on return from ma28a/ad or ma28b/bd, w(1) will give
c     an estimate (an upper bound) of the increase in size of elements
c     encountered during the decomposition. if the matrix is well
c     scaled, then a high value for w(1), relative to the largest entry
c     in the input matrix, indicates that the lu decomposition may be
c     inaccurate and the user should be wary of his results and perhaps
c     increase u for subsequent runs.  we would like to emphasise that
c     this value only relates to the accuracy of our lu decomposition
c     and gives no indication as to the singularity of the matrix or the
c     accuracy of the solution.  this upper bound can be a significant
c     overestimate particularly if the matrix is badly scaled. if an
c     accurate value for the growth is required, lbig (q.v.) should be
c     set to .true.
c eps,rmin are real variables. if, on entry to ma28b/bd, eps is less
c     than one, then rmin will give the smallest ratio of the pivot to
c     the largest element in the corresponding row of the upper
c     triangular factor thus monitoring the stability of successive
c     factorizations. if rmin becomes very large and w(1) from
c     ma28b/bd is also very large, it may be advisable to perform a
c     new decomposition using ma28a/ad.
c resid is a real variable which on exit from ma28c/cd gives the value
c     of the maximum residual over all the equations unsatisfied because
c     of dependency (zero pivots).
c irncp,icncp are integer variables which monitor the adequacy of "elbow
c     room" in irn and a/icn respectively. if either is quite large (say
c     greater than n/10), it will probably pay to increase the size of
c     the corresponding array for subsequent runs. if either is very low
c     or zero then one can perhaps save storage by reducing the size of
c     the corresponding array.
c minirn,minicn are integer variables which, in the event of a
c     successful return (iflag ge 0 or iflag=-14) give the minimum size
c     of irn and a/icn respectively which would enable a successful run
c     on an identical matrix. on an exit with iflag equal to -5, minicn
c     gives the minimum value of icn for success on subsequent runs on
c     an identical matrix. in the event of failure with iflag= -6, -4,
c     -3, -2, or -1, then minicn and minirn give the minimum value of
c     licn and lirn respectively which would be required for a
c     successful decomposition up to the point at which the failure
c     occurred.
c irank is an integer variable which gives an upper bound on the rank of
c     the matrix.
c abort1 is a logical variable with default value .true.  if abort1 is
c     set to .false.  then ma28a/ad will decompose structurally singular
c     matrices (including rectangular ones).
c abort2 is a logical variable with default value .true.  if abort2 is
c     set to .false. then ma28a/ad will decompose numerically singular
c     matrices.
c idisp is an integer array of length 2. on output from ma28a/ad, the
c     indices of the diagonal blocks of the factors lie in positions
c     idisp(1) to idisp(2) of a/icn. this array must be preserved
c     between a call to ma28a/ad and subsequent calls to ma28b/bd,
c     ma28c/cd or ma28i/id.
c tol is a real variable.  if it is set to a positive value, then any
c     non-zero whose modulus is less than tol will be dropped from the
c     factorization.  the factorization will then require less storage
c     but will be inaccurate.  after a run of ma28a/ad with tol positive
c     it is not possible to use ma28b/bd and the user is recommended to
c     use ma28i/id to obtain the solution.  the default value for tol is
c     0.0.
c themax is a real variable.  on exit from ma28a/ad, it will hold the
c     largest entry of the original matrix.
c big is a real variable. if lbig has been set to .true., big will hold
c     the largest entry encountered during the factorization by ma28a/ad
c     or ma28b/bd.
c dxmax is a real variable. on exit from ma28i/id, dxmax will be set to
c     the largest component of the solution.
c errmax is a real variable.  on exit from ma28i/id, if maxit is
c     positive, errmax will be set to the largest component in the
c     estimate of the error.
c dres is a real variable.  on exit from ma28i/id, if maxit is positive,
c     dres will be set to the largest component of the residual.
c cgce is a real variable. it is used by ma28i/id to check the
c     convergence rate.  if the ratio of successive corrections is
c     not less than cgce then we terminate since the convergence
c     rate is adjudged too slow.
c ndrop is an integer variable. if tol has been set positive, on exit
c     from ma28a/ad, ndrop will hold the number of entries dropped from
c     the data structure.
c maxit is an integer variable. it is the maximum number of iterations
c     performed by ma28i/id. it has a default value of 16.
c noiter is an integer variable. it is set by ma28i/id to the number of
c     iterative refinement iterations actually used.
c nsrch is an integer variable. if nsrch is set to a value less than n,
c     then a different pivot option will be employed by ma28a/ad.  this
c     may result in different fill-in and execution time for ma28a/ad.
c     if nsrch is less than or equal to n, the workspace array iw can be
c     reduced in length.  the default value for nsrch is 32768.
c istart is an integer variable. if istart is set to a value other than
c     zero, then the user must supply an estimate of the solution to
c     ma28i/id.  the default value for istart is zero.
c lbig is a logical variable. if lbig is set to .true., the value of the
c     largest element encountered in the factorization by ma28a/ad or
c     ma28b/bd is returned in big.  setting lbig to .true.  will
c     increase the time for ma28a/ad marginally and that for ma28b/bd
c     by about 20%.  the default value for lbig is .false.
c
      DOUBLE PRECISION EPS, RMIN, RESID, TOL, THEMAX, BIG, DXMAX,
     * ERRMAX, DRES, CGCE
      LOGICAL LBLOCK, GROW, ABORT1, ABORT2, LBIG
      COMMON /MA28ED/ LP, MP, LBLOCK, GROW
      COMMON /MA28FD/ EPS, RMIN, RESID, IRNCP, ICNCP, MINIRN, MINICN,
     * IRANK, ABORT1, ABORT2
c     common /ma28gd/ idisp(2)
      COMMON /MA28HD/ TOL, THEMAX, BIG, DXMAX, ERRMAX, DRES, CGCE,
     * NDROP, MAXIT, NOITER, NSRCH, ISTART, LBIG
      DATA EPS /1.0D-4/, TOL /0.0D0/, CGCE /0.5D0/
      DATA MAXIT /16/
      DATA LP /6/, MP /6/, NSRCH /32768/, ISTART /0/
      DATA LBLOCK /.TRUE./, GROW /.TRUE./, LBIG /.FALSE./
      DATA ABORT1 /.TRUE./, ABORT2 /.TRUE./
      END
      BLOCK DATA MABLD2
c although all common block variables do not have default values,
c     we comment on all the common block variables here.
c
c common block ma30e/ed holds control parameters ....
c     common /ma30ed/ lp, abort1, abort2, abort3
c the integer lp is the unit number to which the error messages are
c     sent. lp has a default value of 6.  this default value can be
c     reset by the user, if desired.  a value of 0 suppresses all
c     messages.
c the logical variables abort1,abort2,abort3 are used to control the
c     conditions under which the subroutine will terminate.
c if abort1 is .true. then the subroutine will exit  immediately on
c     detecting structural singularity.
c if abort2 is .true. then the subroutine will exit immediately on
c     detecting numerical singularity.
c if abort3 is .true. then the subroutine will exit immediately when
c     the available space in a/icn is filled up by the previously
c     decomposed, active, and undecomposed parts of the matrix.
c the default values for abort1,abort2,abort3 are set to .true.,.true.
c     and .false. respectively.
c
c the variables in the common block ma30f/fd are used to provide the
c     user with information on the decomposition.
c     common /ma30fd/ irncp, icncp, irank, minirn, minicn
c irncp and icncp are integer variables used to monitor the adequacy
c     of the allocated space in arrays irn and a/icn respectively, by
c     taking account of the number of data management compresses
c     required on these arrays. if irncp or icncp is fairly large (say
c     greater than n/10), it may be advantageous to increase the size
c     of the corresponding array(s).  irncp and icncp are initialized
c     to zero on entry to ma30a/ad and are incremented each time the
c     compressing routine ma30d/dd is entered.
c icncp is the number of compresses on a/icn.
c irncp is the number of compresses on irn.
c irank is an integer variable which gives an estimate (actually an
c     upper bound) of the rank of the matrix. on an exit with iflag
c     equal to 0, this will be equal to n.
c minirn is an integer variable which, after a successful call to
c     ma30a/ad, indicates the minimum length to which irn can be
c     reduced while still permitting a successful decomposition of the
c     same matrix. if, however, the user were to decrease the length
c     of irn to that size, the number of compresses (irncp) may be
c     very high and quite costly. if lirn is not large enough to begin
c     the decomposition on a diagonal block, minirn will be equal to
c     the value required to continue the decomposition and iflag will
c     be set to -3 or -6. a value of lirn slightly greater than this
c     (say about n/2) will usually provide enough space to complete
c     the decomposition on that block. in the event of any other
c     failure minirn gives the minimum size of irn required for a
c     successful decomposition up to that point.
c minicn is an integer variable which after a successful call to
c     ma30a/ad, indicates the minimum size of licn required to enable
c     a successful decomposition. in the event of failure with iflag=
c     -5, minicn will, if abort3 is left set to .false., indicate the
c     minimum length that would be sufficient to prevent this error in
c     a subsequent run on an identical matrix. again the user may
c     prefer to use a value of icn slightly greater than minicn for
c     subsequent runs to avoid too many conpresses (icncp). in the
c     event of failure with iflag equal to any negative value except
c     -4, minicn will give the minimum length to which licn could be
c     reduced to enable a successful decomposition to the point at
c     which failure occurred.  notice that, on a successful entry
c     idisp(2) gives the amount of space in a/icn required for the
c     decomposition while minicn will usually be slightly greater
c     because of the need for "elbow room".  if the user is very
c     unsure how large to make licn, the variable minicn can be used
c     to provide that information. a preliminary run should be
c     performed with abort3 left set to .false. and licn about 3/2
c     times as big as the number of non-zeros in the original matrix.
c     unless the initial problem is very sparse (when the run will be
c     successful) or fills in extremely badly (giving an error return
c     with iflag equal to -4), an error return with iflag equal to -5
c     should result and minicn will give the amount of space required
c     for a successful decomposition.
c
c common block ma30g/gd is used by the ma30b/bd entry only.
c     common /ma30gd/ eps, rmin
c eps is a real/double precision variable. it is used to test for
c     small pivots. its default value is 1.0e-4 (1.0d-4 in d version).
c     if the user sets eps to any value greater than 1.0, then no
c     check is made on the size of the pivots. although the absence of
c     such a check would fail to warn the user of bad instability, its
c     absence will enable ma30b/bd to run slightly faster. an  a
c     posteriori  check on the stability of the factorization can be
c     obtained from mc24a/ad.
c rmin is a real/double precision variable which gives the user some
c     information about the stability of the decomposition.  at each
c     stage of the lu decomposition the magnitude of the pivot apiv
c     is compared with the largest off-diagonal entry currently in its
c     row (row of u), rowmax say. if the ratio
c                       min (apiv/rowmax)
c     where the minimum is taken over all the rows, is less than eps
c     then rmin is set to this minimum value and iflag is returned
c     with the value +i where i is the row in which this minimum
c     occurs.  if the user sets eps greater than one, then this test
c     is not performed. in this case, and when there are no small
c     pivots rmin will be set equal to eps.
c
c common block ma30h/hd is used by ma30c/cd only.
c     common /ma30hd/ resid
c resid is a real/double precision variable. in the case of singular
c     or rectangular matrices its final value will be equal to the
c     maximum residual for the unsatisfied equations; otherwise its
c     value will be set to zero.
c
c common  block ma30i/id controls the use of drop tolerances, the
c     modified pivot option and the the calculation of the largest
c     entry in the factorization process. this common block was added
c     to the ma30 package in february, 1983.
c     common /ma30id/ tol, big, ndrop, nsrch, lbig
c tol is a real/double precision variable.  if it is set to a positive
c     value, then ma30a/ad will drop from the factors any non-zero
c     whose modulus is less than tol.  the factorization will then
c     require less storage but will be inaccurate.  after a run of
c     ma30a/ad where entries have been dropped, ma30b/bd  should not
c     be called.  the default value for tol is 0.0.
c big is a real/double precision variable.  if lbig has been set to
c     .true., big will be set to the largest entry encountered during
c     the factorization.
c ndrop is an integer variable. if tol has been set positive, on exit
c     from ma30a/ad, ndrop will hold the number of entries dropped
c     from the data structure.
c nsrch is an integer variable. if nsrch is set to a value less than
c     or equal to n, then a different pivot option will be employed by
c     ma30a/ad.  this may result in different fill-in and execution
c     time for ma30a/ad. if nsrch is less than or equal to n, the
c     workspace arrays lastc and nextc are not referenced by ma30a/ad.
c     the default value for nsrch is 32768.
c lbig is a logical variable. if lbig is set to .true., the value of
c     the largest entry encountered in the factorization by ma30a/ad
c     is returned in big.  setting lbig to .true.  will marginally
c     increase the factorization time for ma30a/ad and will increase
c     that for ma30b/bd by about 20%.  the default value for lbig is
c     .false.
c
      DOUBLE PRECISION EPS, RMIN, TOL, BIG
      LOGICAL ABORT1, ABORT2, ABORT3, LBIG
      COMMON /MA30ED/ LP, ABORT1, ABORT2, ABORT3
      COMMON /MA30GD/ EPS, RMIN
      COMMON /MA30ID/ TOL, BIG, NDROP, NSRCH, LBIG
      DATA EPS /1.0D-4/, TOL /0.0D0/, BIG /0.0D0/
      DATA LP /6/, NSRCH /32768/
      DATA LBIG /.FALSE./
      DATA ABORT1 /.TRUE./, ABORT2 /.TRUE./, ABORT3 /.FALSE./
      END
      BLOCK DATA MABLD3
      LOGICAL ABORT
      COMMON /MC23BD/ LP,NUMNZ,NUM,LARGE,ABORT
      DATA LP/6/,ABORT/.FALSE./
      END
