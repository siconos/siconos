c-----------------------------------------------------------------------
c Demonstration program for DLSODPK.
c ODE system from ns-species interaction pde in 2 dimensions.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c-----------------------------------------------------------------------
c This program solves a stiff ODE system that arises from a system
c of partial differential equations.  The PDE system is a food web
c population model, with predator-prey interaction and diffusion on
c the unit square in two dimensions.  The dependent variable vector is
c
c         1   2        ns
c   c = (c , c , ..., c  )
c
c and the PDEs are as follows:
c
c     i               i      i
c   dc /dt  =  d(i)*(c    + c   )  +  f (x,y,c)  (i=1,...,ns)
c                     xx     yy        i
c
c where
c                  i          ns         j
c   f (x,y,c)  =  c *(b(i) + sum a(i,j)*c )
c    i                       j=1
c
c The number of species is ns = 2*np, with the first np being prey and
c the last np being predators.  The coefficients a(i,j), b(i), d(i) are:
c
c   a(i,i) = -a  (all i)
c   a(i,j) = -g  (i .le. np, j .gt. np)
c   a(i,j) =  e  (i .gt. np, j .le. np)
c   b(i) =  b*(1 + alpha*x*y)  (i .le. np)
c   b(i) = -b*(1 + alpha*x*y)  (i .gt. np)
c   d(i) = dprey  (i .le. np)
c   d(i) = dpred  (i .gt. np)
c
c The various scalar parameters are set in subroutine setpar.
c
c The boundary conditions are: normal derivative = 0.
c A polynomial in x and y is used to set the initial conditions.
c
c The PDEs are discretized by central differencing on a mx by my mesh.
c
c The ODE system is solved by DLSODPK using method flag values
c mf = 10, 21, 22, 23, 24, 29.  The final time is tmax = 10, except
c that for mf = 10 it is tmax = 1.0d-3 because the problem is stiff,
c and for mf = 23 and 24 it is tmax = 2 because the lack of symmetry
c in the problem makes these methods more costly.
c
c Two preconditioner matrices are used.  One uses a fixed number of
c Gauss-Seidel iterations based on the diffusion terms only.
c The other preconditioner is a block-diagonal matrix based on
c the partial derivatives of the interaction terms f only, using
c block-grouping (computing only a subset of the ns by ns blocks).
c For mf = 21 and 22, these two preconditioners are applied on
c the left and right, respectively, and for mf = 23 and 24 the product
c of the two is used as the one preconditioner matrix.
c For mf = 29, the inverse of the product is applied.
c
c Two output files are written: one with the problem description and
c and performance statistics on unit 6, and one with solution profiles
c at selected output times (for mf = 22 only) on unit 8.
c-----------------------------------------------------------------------
c Note: In addition to the main program and 10 subroutines
c given below, this program requires the LINPACK subroutines
c DGEFA and DGESL, and the BLAS routine DAXPY.
c-----------------------------------------------------------------------
c Reference:
c     Peter N. Brown and Alan C. Hindmarsh,
c     Reduced Storage Matrix Methods in Stiff ODE Systems,
c     J. Appl. Math. & Comp., 31 (1989), pp. 40-91;
c     Also LLNL Report UCRL-95088, Rev. 1, June 1987.
c-----------------------------------------------------------------------
      external fweb, jacbg, solsbg
      integer ns, mx, my, mxns,
     1        mp, mq, mpsq, itmax,
     2        meshx,meshy,ngx,ngy,ngrp,mxmp,jgx,jgy,jigx,jigy,jxr,jyr
      integer i, imod3, iopt, iout, istate, itask, itol, iwork,
     1   jacflg, jpre, leniw, lenrw, liw, lrw, mf,
     2   ncfl, ncfn, neq, nfe, nfldif, nfndif, nli, nlidif, nni, nnidif,
     3   nout, npe, nps, nqu, nsdif, nst
      double precision aa, ee, gg, bb, dprey, dpred,
     1     ax, ay, acoef, bcoef, dx, dy, alph, diff, cox, coy,
     2     uround, srur
      double precision avdim, atol, cc, hu, rcfl, rcfn, dumach,
     1   rtol, rwork, t, tout
c
c The problem Common blocks below allow for up to 20 species,
c up to a 50x50 mesh, and up to a 20x20 group structure.
      common /pcom0/ aa, ee, gg, bb, dprey, dpred
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
      common /pcom2/ uround, srur, mp, mq, mpsq, itmax
      common /pcom3/ meshx, meshy, ngx, ngy, ngrp, mxmp,
     2   jgx(21), jgy(21), jigx(50), jigy(50), jxr(20), jyr(20)
c
c The dimension of cc below must be .ge. 2*neq, where neq = ns*mx*my.
c The dimension lrw of rwork must be .ge. 17*neq + ns*ns*ngrp + 61,
c and the dimension liw of iwork must be .ge. 35 + ns*ngrp.
      dimension cc(576), rwork(5213), iwork(67)
      data lrw/5213/, liw/67/
c
      open (unit=6, file='demout', status='new')
      open (unit=8, file='ccout', status='new')
c
      ax = 1.0d0
      ay = 1.0d0
c
c Call setpar to set problem parameters.
      call setpar
c
c Set remaining problem parameters.
      neq = ns*mx*my
      mxns = mx*ns
      dx = ax/(mx-1)
      dy = ay/(my-1)
      do 10 i = 1,ns
        cox(i) = diff(i)/dx**2
 10     coy(i) = diff(i)/dy**2
c
c Write heading.
      write(6,20)ns, mx,my,neq
 20   format(' Demonstration program for DLSODPK package'//
     1   ' Food web problem with ns species, ns =',i4/
     2   ' Predator-prey interaction and diffusion on a 2-d square'//
     3   ' Mesh dimensions (mx,my) =',2i4/
     4   ' Total system size is neq =',i7//)
      write(6,25) aa,ee,gg,bb,dprey,dpred,alph
 25   format(' Matrix parameters:  a =',d12.4,'   e =',d12.4,
     1   '   g =',d12.4/20x,' b =',d12.4//
     2   ' Diffusion coefficients: dprey =',d12.4,'   dpred =',d12.4/
     3   ' Rate parameter alpha =',d12.4//)
c
c Set remaining method parameters.
      jpre = 3
      jacflg = 1
      iwork(3) = jpre
      iwork(4) = jacflg
      iopt = 0
      mp = ns
      mq = mx*my
      mpsq = ns*ns
      uround = dumach()
      srur = sqrt(uround)
      meshx = mx
      meshy = my
      mxmp = meshx*mp
      ngx = 2
      ngy = 2
      ngrp = ngx*ngy
      call gset (meshx, ngx, jgx, jigx, jxr)
      call gset (meshy, ngy, jgy, jigy, jyr)
      iwork(1) = mpsq*ngrp
      iwork(2) = mp*ngrp
      itmax = 5
      itol = 1
      rtol = 1.0d-5
      atol = rtol
      itask = 1
      write(6,30)ngrp,ngx,ngy,itmax,rtol,atol
 30   format(' Preconditioning uses interaction-only block-diagonal',
     1   ' matrix'/' with block-grouping, and Gauss-Seidel iterations'//
     2   ' Number of diagonal block groups = ngrp =',i4,
     3   '   (ngx by ngy, ngx =',i2,'  ngy =',i2,' )'//
     4   ' G-S preconditioner uses itmax iterations, itmax =',i3//
     5   ' Tolerance parameters: rtol =',d10.2,'   atol =',d10.2)
c
c
c Loop over mf values 10, 21, ..., 24, 29.
c
      do 90 mf = 10,29
      if (mf .gt. 10 .and. mf .lt. 21) go to 90
      if (mf .gt. 24 .and. mf .lt. 29) go to 90
      write(6,40)mf
 40   format(//80('-')//' Solution with mf =',i3//
     1   '   t       nstep  nfe  nni  nli  npe  nq',
     2   4x,'h          avdim    ncf rate    lcf rate')
c
      t = 0.0d0
      tout = 1.0d-8
      nout = 18
      if (mf .eq. 10) nout = 6
      if (mf .eq. 23 .or. mf .eq. 24) nout = 10
      call cinit (cc)
      if (mf .eq. 22) call outweb (t, cc, ns, mx, my, 8)
      istate = 1
      nli = 0
      nni = 0
      ncfn = 0
      ncfl = 0
      nst = 0
c
c Loop over output times and call DLSODPK.
c
      do 70 iout = 1,nout
        call dlsodpk (fweb, neq, cc, t, tout, itol, rtol, atol, itask,
     1         istate, iopt, rwork, lrw, iwork, liw, jacbg, solsbg, mf)
        nsdif = iwork(11) - nst
        nst = iwork(11)
        nfe = iwork(12)
        npe = iwork(13)
        nnidif = iwork(19) - nni
        nni = iwork(19)
        nlidif = iwork(20) - nli
        nli = iwork(20)
        nfndif = iwork(22) - ncfn
        ncfn = iwork(22)
        nfldif = iwork(23) - ncfl
        ncfl = iwork(23)
        nqu = iwork(14)
        hu = rwork(11)
        avdim = 0.0d0
        rcfn = 0.0d0
        rcfl = 0.0d0
        if (nnidif .gt. 0) avdim = real(nlidif)/real(nnidif)
        if (nsdif .gt. 0) rcfn = real(nfndif)/real(nsdif)
        if (nnidif .gt. 0) rcfl = real(nfldif)/real(nnidif)
        write(6,50)t,nst,nfe,nni,nli,npe,nqu,hu,avdim,rcfn,rcfl
 50     format(d10.2,i5,i6,3i5,i4,2d11.2,d10.2,d12.2)
        imod3 = iout - 3*(iout/3)
        if (mf .eq. 22 .and. imod3 .eq. 0) call outweb (t,cc,ns,mx,my,8)
        if (istate .eq. 2) go to 65
        write(6,60)t
 60     format(//' final time reached =',d12.4//)
        go to 75
 65     continue
        if (tout .gt. 0.9d0) tout = tout + 1.0d0
        if (tout .lt. 0.9d0) tout = tout*10.0d0
 70     continue
c
 75   continue
      nst = iwork(11)
      nfe = iwork(12)
      npe = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nni = iwork(19)
      nli = iwork(20)
      nps = iwork(21)
      if (nni .gt. 0) avdim = real(nli)/real(nni)
      ncfn = iwork(22)
      ncfl = iwork(23)
      write (6,80) lenrw,leniw,nst,nfe,npe,nps,nni,nli,avdim,
     1               ncfn,ncfl
 80   format(//' Final statistics for this run:'/
     1   ' rwork size =',i8,'   iwork size =',i6/
     2   ' number of time steps            =',i5/
     3   ' number of f evaluations         =',i5/
     4   ' number of preconditioner evals. =',i5/
     4   ' number of preconditioner solves =',i5/
     5   ' number of nonlinear iterations  =',i5/
     5   ' number of linear iterations     =',i5/
     6   ' average subspace dimension  =',f8.4/
     7   i5,' nonlinear conv. failures,',i5,' linear conv. failures')
c
 90   continue
      stop
c------  end of main program for DLSODPK demonstration program ----------
      end

      subroutine setpar
c-----------------------------------------------------------------------
c This routine sets the problem parameters.
c It set ns, mx, my, and problem coefficients acoef, bcoef, diff, alph,
c using parameters np, aa, ee, gg, bb, dprey, dpred.
c-----------------------------------------------------------------------
      integer ns, mx, my, mxns
      integer i, j, np
      double precision aa, ee, gg, bb, dprey, dpred,
     1     ax, ay, acoef, bcoef, dx, dy, alph, diff, cox, coy
      common /pcom0/ aa, ee, gg, bb, dprey, dpred
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
c
      np = 3
      mx = 6
      my = 6
      aa = 1.0d0
      ee = 1.0d4
      gg = 0.5d-6
      bb = 1.0d0
      dprey = 1.0d0
      dpred = 0.5d0
      alph = 1.0d0
      ns = 2*np
      do 70 j = 1,np
        do 60 i = 1,np
          acoef(np+i,j) = ee
          acoef(i,np+j) = -gg
 60       continue
        acoef(j,j) = -aa
        acoef(np+j,np+j) = -aa
        bcoef(j) = bb
        bcoef(np+j) = -bb
        diff(j) = dprey
        diff(np+j) = dpred
 70     continue
c
      return
c------------  end of subroutine setpar  -------------------------------
      end

      subroutine gset (m, ng, jg, jig, jr)
c-----------------------------------------------------------------------
c This routine sets arrays jg, jig, and jr describing
c a uniform partition of (1,2,...,m) into ng groups.
c-----------------------------------------------------------------------
      integer m, ng, jg, jig, jr
      dimension jg(*), jig(*), jr(*)
      integer ig, j, len1, mper, ngm1
c
      mper = m/ng
      do 10 ig = 1,ng
 10     jg(ig) = 1 + (ig - 1)*mper
      jg(ng+1) = m + 1
c
      ngm1 = ng - 1
      len1 = ngm1*mper
      do 20 j = 1,len1
 20     jig(j) = 1 + (j-1)/mper
      len1 = len1 + 1
      do 25 j = len1,m
 25     jig(j) = ng
c
      do 30 ig = 1,ngm1
 30     jr(ig) = 0.5d0 + (ig - 0.5d0)*mper
      jr(ng) = 0.5d0*(1 + ngm1*mper + m)
c
      return
c------------  end of subroutine gset  ---------------------------------
      end

      subroutine cinit (cc)
c-----------------------------------------------------------------------
c This routine computes and loads the vector of initial values.
c-----------------------------------------------------------------------
      double precision cc
      dimension cc(*)
      integer ns, mx, my, mxns
      integer i, ici, ioff, iyoff, jx, jy
      double precision ax,ay,acoef,bcoef,dx,dy,alph,diff,cox,coy
      double precision argx, argy, x, y
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
c
        do 20 jy = 1,my
          y = (jy-1)*dy
          argy = 16.0d0*y*y*(ay-y)*(ay-y)
          iyoff = mxns*(jy-1)
          do 10 jx = 1,mx
            x = (jx-1)*dx
            argx = 16.0d0*x*x*(ax-x)*(ax-x)
            ioff = iyoff + ns*(jx-1)
            do 5 i = 1,ns
              ici = ioff + i
              cc(ici) = 10.0d0 + i*argx*argy
  5           continue
 10         continue
 20       continue
      return
c------------  end of subroutine cinit  --------------------------------
      end

      subroutine outweb (t, c, ns, mx, my, lun)
c-----------------------------------------------------------------------
c This routine prints the values of the individual species densities
c at the current time t.  The write statements use unit lun.
c-----------------------------------------------------------------------
      integer ns, mx, my, lun
      double precision t, c
      dimension c(ns,mx,my)
      integer i, jx, jy
c
      write(lun,10) t
 10   format(/80('-')/30x,'At time t = ',d16.8/80('-') )
c
      do 40 i = 1,ns
        write(lun,20) i
 20     format(' the species c(',i2,') values are:')
        do 30 jy = my,1,-1
          write(lun,25) (c(i,jx,jy),jx=1,mx)
 25       format(6(1x,g12.6))
 30       continue
        write(lun,35)
 35     format(80('-'),/)
 40     continue
c
      return
c------------  end of subroutine outweb  -------------------------------
      end

      subroutine fweb (neq, t, cc, cdot)
c-----------------------------------------------------------------------
c This routine computes the derivative of cc and returns it in cdot.
c The interaction rates are computed by calls to webr, and these are
c saved in cc(neq+1),...,cc(2*neq) for use in preconditioning.
c-----------------------------------------------------------------------
      integer neq
      double precision t, cc, cdot
      dimension cc(neq), cdot(neq)
      integer ns, mx, my, mxns
      integer i, ic, ici, idxl, idxu, idyl, idyu, iyoff, jx, jy
      double precision ax,ay,acoef,bcoef,dx,dy,alph,diff,cox,coy
      double precision dcxli, dcxui, dcyli, dcyui, x, y
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
c
      do 100 jy = 1,my
        y = (jy-1)*dy
        iyoff = mxns*(jy-1)
        idyu = mxns
        if (jy .eq. my) idyu = -mxns
        idyl = mxns
        if (jy .eq. 1) idyl = -mxns
        do 90 jx = 1,mx
          x = (jx-1)*dx
          ic = iyoff + ns*(jx-1) + 1
c Get interaction rates at one point (x,y).
          call webr (x, y, t, cc(ic), cc(neq+ic))
          idxu = ns
          if (jx .eq. mx) idxu = -ns
          idxl = ns
          if (jx .eq. 1) idxl = -ns
          do 80 i = 1,ns
            ici = ic + i - 1
c Do differencing in y.
            dcyli = cc(ici) - cc(ici-idyl)
            dcyui = cc(ici+idyu) - cc(ici)
c Do differencing in x.
            dcxli = cc(ici) - cc(ici-idxl)
            dcxui = cc(ici+idxu) - cc(ici)
c Collect terms and load cdot elements.
            cdot(ici) = coy(i)*(dcyui - dcyli) + cox(i)*(dcxui - dcxli)
     1                  + cc(neq+ici)
 80         continue
 90       continue
 100    continue
      return
c------------  end of subroutine fweb  ---------------------------------
      end

      subroutine webr (x, y, t, c, rate)
c-----------------------------------------------------------------------
c This routine computes the interaction rates for the species
c c(1),...,c(ns), at one spatial point and at time t.
c-----------------------------------------------------------------------
      double precision x, y, t, c, rate
      dimension c(*), rate(*)
      integer ns, mx, my, mxns
      integer i
      double precision ax,ay,acoef,bcoef,dx,dy,alph,diff,cox,coy
      double precision fac
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
c
      do 10 i = 1,ns
 10     rate(i) = 0.0d0
      do 15 i = 1,ns
        call daxpy (ns, c(i), acoef(1,i), 1, rate, 1)
 15     continue
      fac = 1.0d0 + alph*x*y
      do 20 i = 1,ns
 20     rate(i) = c(i)*(bcoef(i)*fac + rate(i))
      return
c------------  end of subroutine webr  ---------------------------------
      end

      subroutine jacbg (f, neq, t, cc, ccsv, rewt, f0, f1, hl0,
     1                  bd, ipbd, ier)
c-----------------------------------------------------------------------
c This routine generates part of the block-diagonal part of the
c Jacobian, multiplies by -hl0, adds the identity matrix,
c and calls DGEFA to do LU decomposition of each diagonal block.
c The computation of the diagonal blocks uses the block and grouping
c information in /pcom1/ and /pcom2/.  One block per group is computed.
c The Jacobian elements are generated by difference quotients
c using calls to the routine fbg.
c-----------------------------------------------------------------------
c The two Common blocks below are used for internal communication.
c The variables used are:
c   mp     = size of blocks in block-diagonal preconditioning matrix.
c   mq     = number of blocks in each direction (neq = mp*mq).
c   mpsq   = mp*mp.
c   uround = unit roundoff, generated by a call uround = dumach().
c   srur   = sqrt(uround).
c   meshx  = x mesh size
c   meshy  = y mesh size (mesh is meshx by meshy)
c   ngx    = no. groups in x direction in block-grouping scheme.
c   ngy    = no. groups in y direction in block-grouping scheme.
c   ngrp   = total number of groups = ngx*ngy.
c   mxmp   = meshx*mp.
c   jgx    = length ngx+1 array of group boundaries in x direction.
c            group igx has x indices jx = jgx(igx),...,jgx(igx+1)-1.
c   jigx   = length meshx array of x group indices vs x node index.
c            x node index jx is in x group jigx(jx).
c   jxr    = length ngx array of x indices representing the x groups.
c            the index for x group igx is jx = jxr(igx).
c   jgy, jigy, jyr = analogous arrays for grouping in y direction.
c-----------------------------------------------------------------------
      external f
      integer neq, ipbd, ier
      double precision t, cc, ccsv, rewt, f0, f1, hl0, bd
      dimension cc(neq), ccsv(neq), rewt(neq), f0(neq), f1(neq),
     1          bd(*), ipbd(*)
      integer mp, mq, mpsq, itmax,
     2        meshx,meshy,ngx,ngy,ngrp,mxmp,jgx,jgy,jigx,jigy,jxr,jyr
      integer i, ibd, idiag, if0, if00, ig, igx, igy, iip,
     1   j, jj, jx, jy, n
      double precision uround, srur
      double precision fac, r, r0, dvnorm
c
      common /pcom2/ uround, srur, mp, mq, mpsq, itmax
      common /pcom3/ meshx, meshy, ngx, ngy, ngrp, mxmp,
     1   jgx(21), jgy(21), jigx(50), jigy(50), jxr(20), jyr(20)
c
      n = neq
c
c-----------------------------------------------------------------------
c Make mp calls to fbg to approximate each diagonal block of Jacobian.
c Here cc(neq+1),...,cc(2*neq) contains the base fb value.
c r0 is a minimum increment factor for the difference quotient.
c-----------------------------------------------------------------------
 200  fac = dvnorm (n, f0, rewt)
      r0 = 1000.0d0*abs(hl0)*uround*n*fac
      if (r0 .eq. 0.0d0) r0 = 1.0d0
      ibd = 0
      do 240 igy = 1,ngy
        jy = jyr(igy)
        if00 = (jy - 1)*mxmp
        do 230 igx = 1,ngx
          jx = jxr(igx)
          if0 = if00 + (jx - 1)*mp
          do 220 j = 1,mp
            jj = if0 + j
            r = max(srur*abs(cc(jj)),r0/rewt(jj))
            cc(jj) = cc(jj) + r
            fac = -hl0/r
            call fbg (neq, t, cc, jx, jy, f1)
            do 210 i = 1,mp
 210          bd(ibd+i) = (f1(i) - cc(neq+if0+i))*fac
            cc(jj) = ccsv(jj)
            ibd = ibd + mp
 220        continue
 230      continue
 240    continue
c
c Add identity matrix and do LU decompositions on blocks. --------------
 260  continue
      ibd = 1
      iip = 1
      do 280 ig = 1,ngrp
        idiag = ibd
        do 270 i = 1,mp
          bd(idiag) = bd(idiag) + 1.0d0
 270      idiag = idiag + (mp + 1)
        call dgefa (bd(ibd), mp, mp, ipbd(iip), ier)
        if (ier .ne. 0) go to 290
        ibd = ibd + mpsq
        iip = iip + mp
 280    continue
 290  return
c------------  end of subroutine jacbg  --------------------------------
      end

      subroutine fbg (neq, t, cc, jx, jy, cdot)
c-----------------------------------------------------------------------
c This routine computes one block of the interaction terms of the
c system, namely block (jx,jy), for use in preconditioning.
c-----------------------------------------------------------------------
      integer neq, jx, jy
      double precision t, cc, cdot
      dimension cc(neq), cdot(neq)
      integer ns, mx, my, mxns
      integer iblok, ic
      double precision ax,ay,acoef,bcoef,dx,dy,alph,diff,cox,coy
      double precision x, y
c
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
c
      iblok = jx + (jy-1)*mx
        y = (jy-1)*dy
          x = (jx-1)*dx
          ic = ns*(iblok-1) + 1
          call webr (x, y, t, cc(ic), cdot)
      return
c------------  end of subroutine fbg  ----------------------------------
      end

      subroutine solsbg (n, t, cc, f0, wk, hl0, bd, ipbd, v, lr, ier)
c-----------------------------------------------------------------------
c This routine applies one or two inverse preconditioner matrices
c to the array v, using the interaction-only block-diagonal Jacobian
c with block-grouping, and Gauss-Seidel applied to the diffusion terms.
c When lr = 1 or 3, it calls gs for a Gauss-Seidel approximation
c to ((I-hl0*Jd)-inverse)*v, and stores the result in v.
c When lr = 2 or 3, it computes ((I-hl0*dg/dc)-inverse)*v, using LU
c factors of the blocks in bd, and pivot information in ipbd.
c In both cases, the array v is overwritten with the solution.
c-----------------------------------------------------------------------
      integer n, ipbd, lr, ier
      double precision t, cc, f0, wk, hl0, bd, v
      dimension cc(n), f0(n), wk(n), bd(*), ipbd(*), v(n)
      integer mp, mq, mpsq, itmax,
     2        meshx,meshy,ngx,ngy,ngrp,mxmp,jgx,jgy,jigx,jigy,jxr,jyr
      integer ibd, ig0, igm1, igx, igy, iip, iv, jx, jy
      double precision uround, srur
c
      common /pcom2/ uround, srur, mp, mq, mpsq, itmax
      common /pcom3/ meshx, meshy, ngx, ngy, ngrp, mxmp,
     1   jgx(21), jgy(21), jigx(50), jigy(50), jxr(20), jyr(20)
c
      ier = 0
c
      if (lr.eq.0 .or. lr.eq.1 .or. lr.eq.3) call gs (n, hl0, v, wk)
      if (lr.eq.0 .or. lr.eq.2 .or. lr.eq.3) then
        iv = 1
        do 20 jy = 1,meshy
          igy = jigy(jy)
          ig0 = (igy - 1)*ngx
          do 10 jx = 1,meshx
            igx = jigx(jx)
            igm1 = igx - 1 + ig0
            ibd = 1 + igm1*mpsq
            iip = 1 + igm1*mp
            call dgesl (bd(ibd), mp, mp, ipbd(iip), v(iv), 0)
            iv = iv + mp
 10         continue
 20       continue
        endif
c
      return
c------------  end of subroutine solsbg  -------------------------------
      end

      subroutine gs (n, hl0, z, x)
c-----------------------------------------------------------------------
c This routine performs itmax Gauss-Seidel iterations
c to compute an approximation to P-inverse*z, where P = I - hl0*Jd, and
c Jd represents the diffusion contributions to the Jacobian.
c z contains the answer on return.
c The dimensions below assume ns .le. 20.
c-----------------------------------------------------------------------
      integer n
      double precision hl0, z, x
      dimension z(n), x(n)
      integer ns, mx, my, mxns,
     1        mp, mq, mpsq, itmax
      integer i, ic, ici, iter, iyoff, jx, jy
      double precision ax,ay,acoef,bcoef,dx,dy,alph,diff,cox,coy,
     2     uround, srur
      double precision beta,beta2,cof1,elamda,gamma,gamma2
      dimension beta(20), gamma(20), beta2(20), gamma2(20), cof1(20)
      common /pcom1/ ax, ay, acoef(20,20), bcoef(20), dx, dy, alph,
     1               diff(20), cox(20), coy(20), ns, mx, my, mxns
      common /pcom2/ uround, srur, mp, mq, mpsq, itmax
c
c-----------------------------------------------------------------------
c Write matrix as P = D - L - U.
c Load local arrays beta, beta2, gamma, gamma2, and cof1.
c-----------------------------------------------------------------------
      do 10 i = 1,ns
        elamda = 1.d0/(1.d0 + 2.d0*hl0*(cox(i) + coy(i)))
        beta(i) = hl0*cox(i)*elamda
        beta2(i) = 2.d0*beta(i)
        gamma(i) = hl0*coy(i)*elamda
        gamma2(i) = 2.d0*gamma(i)
        cof1(i) = elamda
 10     continue
c-----------------------------------------------------------------------
c Begin iteration loop.
c Load array x with (D-inverse)*z for first iteration.
c-----------------------------------------------------------------------
      iter = 1
c
      do 50 jy = 1,my
        iyoff = mxns*(jy-1)
        do 40 jx = 1,mx
          ic = iyoff + ns*(jx-1)
          do 30 i = 1,ns
            ici = ic + i
            x(ici) = cof1(i)*z(ici)
            z(ici) = 0.d0
 30         continue
 40       continue
 50     continue
      go to 160
c-----------------------------------------------------------------------
c Calculate (D-inverse)*U*x.
c-----------------------------------------------------------------------
 70   continue
      iter = iter + 1
      jy = 1
        jx = 1
        ic = ns*(jx-1)
        do 75 i = 1,ns
          ici = ic + i
 75       x(ici) = beta2(i)*x(ici+ns) + gamma2(i)*x(ici+mxns)
        do 85 jx = 2,mx-1
          ic = ns*(jx-1)
          do 80 i = 1,ns
            ici = ic + i
 80         x(ici) = beta(i)*x(ici+ns) + gamma2(i)*x(ici+mxns)
 85       continue
        jx = mx
        ic = ns*(jx-1)
        do 90 i = 1,ns
          ici = ic + i
 90       x(ici) = gamma2(i)*x(ici+mxns)
      do 115 jy = 2,my-1
        iyoff = mxns*(jy-1)
          jx = 1
          ic = iyoff
          do 95 i = 1,ns
            ici = ic + i
 95         x(ici) = beta2(i)*x(ici+ns) + gamma(i)*x(ici+mxns)
          do 105 jx = 2,mx-1
            ic = iyoff + ns*(jx-1)
            do 100 i = 1,ns
              ici = ic + i
 100          x(ici) = beta(i)*x(ici+ns) + gamma(i)*x(ici+mxns)
 105        continue
          jx = mx
          ic = iyoff + ns*(jx-1)
          do 110 i = 1,ns
            ici = ic + i
 110        x(ici) = gamma(i)*x(ici+mxns)
 115      continue
      jy = my
      iyoff = mxns*(jy-1)
        jx = 1
        ic = iyoff
        do 120 i = 1,ns
          ici = ic + i
 120      x(ici) = beta2(i)*x(ici+ns)
        do 130 jx = 2,mx-1
          ic = iyoff + ns*(jx-1)
          do 125 i = 1,ns
            ici = ic + i
 125      x(ici) = beta(i)*x(ici+ns)
 130      continue
        jx = mx
        ic = iyoff + ns*(jx-1)
        do 135 i = 1,ns
          ici = ic + i
 135      x(ici) = 0.0d0
c-----------------------------------------------------------------------
c Calculate (I - (D-inverse)*L)-inverse * x.
c-----------------------------------------------------------------------
 160  continue
      jy = 1
        do 175 jx = 2,mx-1
          ic = ns*(jx-1)
          do 170 i = 1,ns
            ici = ic + i
 170        x(ici) = x(ici) + beta(i)*x(ici-ns)
 175      continue
        jx = mx
        ic = ns*(jx-1)
        do 180 i = 1,ns
          ici = ic + i
 180      x(ici) = x(ici) + beta2(i)*x(ici-ns)
      do 210 jy = 2,my-1
        iyoff = mxns*(jy-1)
          jx = 1
          ic = iyoff
          do 185 i = 1,ns
            ici = ic + i
 185        x(ici) = x(ici) + gamma(i)*x(ici-mxns)
          do 200 jx = 2,mx-1
            ic = iyoff + ns*(jx-1)
            do 195 i = 1,ns
              ici = ic + i
              x(ici) = (x(ici) + beta(i)*x(ici-ns))
     1             + gamma(i)*x(ici-mxns)
 195          continue
 200        continue
            jx = mx
            ic = iyoff + ns*(jx-1)
            do 205 i = 1,ns
              ici = ic + i
              x(ici) = (x(ici) + beta2(i)*x(ici-ns))
     1             + gamma(i)*x(ici-mxns)
 205          continue
 210        continue
      jy = my
      iyoff = mxns*(jy-1)
        jx = 1
        ic = iyoff
        do 215 i = 1,ns
          ici = ic + i
 215      x(ici) = x(ici) + gamma2(i)*x(ici-mxns)
        do 225 jx = 2,mx-1
          ic = iyoff + ns*(jx-1)
          do 220 i = 1,ns
            ici = ic + i
            x(ici) = (x(ici) + beta(i)*x(ici-ns))
     1           + gamma2(i)*x(ici-mxns)
 220        continue
 225      continue
        jx = mx
        ic = iyoff + ns*(jx-1)
        do 230 i = 1,ns
          ici = ic + i
          x(ici) = (x(ici) + beta2(i)*x(ici-ns))
     1         + gamma2(i)*x(ici-mxns)
 230      continue
c-----------------------------------------------------------------------
c Add increment x to z.
c-----------------------------------------------------------------------
      do 300 i = 1,n
 300    z(i) = z(i) + x(i)
c
      if (iter .lt. itmax) go to 70
      return
c------------  end of subroutine gs  -----------------------------------
      end
