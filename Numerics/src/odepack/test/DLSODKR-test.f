c-----------------------------------------------------------------------
c Demonstration program for the DLSODKR package.
c This is the version of 15 June 2001.
c
c This version is in double precision.
c
c An ODE system is generated from the following 2-species diurnal
c kinetics advection-diffusion PDE system in 2 space dimensions:
c
c dc(i)/dt = Kh*(d/dx)**2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)
c                 + Ri(c1,c2,t)      for i = 1,2,   where
c   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
c   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
c   Kv(z) = Kv0*exp(z/5) ,
c Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
c vary diurnally.  The species are oxygen singlet and ozone.
c The problem is posed on the square
c   0 .le. x .le. 20,    30 .le. z .le. 50   (all in km),
c with homogeneous Neumann boundary conditions, and for time t in
c   0 .le. t .le. 86400 sec (1 day).
c
c The PDE system is treated by central differences on a uniform
c 10 x 10 mesh, with simple polynomial initial profiles.
c
c The problem is solved with DLSODKR, with the BDF/GMRES method and
c the block-diagonal part of the Jacobian as a left preconditioner.
c At intervals of 7200 sec (2 hrs), output includes sample values
c of c1, c2, and c2tot = total of all c2 values.
c
c Roots of the function g = d(c2tot)/dt are found, i.e. the points at
c which the total c2 (ozone) is stationary.
c
c Note: The preconditioner routines call LINPACK routines DGEFA, DGESL,
c and the BLAS routines DCOPY and DSCAL.
c-----------------------------------------------------------------------
      external fdem, jacbd, solbd, gdem
      integer mx, mz, mm, iout, istate, iwork, jacflg, jpre, jroot,
     1   jx, jz, leniw, lenrw, liw, lrw, mf, neq, nst, nsfi, nfe,
     2   nge, npe, njev, nps, nni, nli, ncfn, ncfl   
      double precision q1,q2,q3,q4, a3,a4, om, c3, dz, hdco, vdco, haco,
     1   dkh, vel, dkv0, halfda, pi, twohr, rtol, floor,
     2   dx, atol, t, tout, x, cx, z, cz, y, rwork, c2tot, avdim
      common /pcom/ q1,q2,q3,q4,a3,a4,om,c3,dz,hdco,vdco,haco,mx,mz,mm
      dimension y(2,10,10), rwork(4264), iwork(235)
      data dkh/4.0d-6/, vel/0.001d0/, dkv0/1.0d-8/, halfda/4.32d4/,
     1  pi/3.1415926535898d0/, twohr/7200.0d0/, rtol/1.0d-5/,
     2  floor/100.0d0/, lrw/4264/, liw/235/, mf/22/, jpre/1/, jacflg/1/
c
c Load Common block of problem parameters.
      mx = 10
      mz = 10
      mm = mx*mz
      q1 = 1.63d-16
      q2 = 4.66d-16
      a3 = 22.62d0
      a4 = 7.601d0
      om = pi/halfda
      c3 = 3.7d16
      dx = 20.0d0/(mx - 1.0d0)
      dz = 20.0d0/(mz - 1.0d0)
      hdco = dkh/dx**2
      haco = vel/(2.0d0*dx)
      vdco = (1.0d0/dz**2)*dkv0
c Set other input arguments.
      atol = rtol*floor
      neq = 2*mx*mz
      iwork(1) = 8*mx*mz
      iwork(2) = neq
      iwork(3) = jpre
      iwork(4) = jacflg
      t = 0.0d0
      tout = 0.0d0
      istate = 1
c Set initial profiles.
      do 20 jz = 1,mz
        z = 30.0d0 + (jz - 1.0d0)*dz
        cz = (0.1d0*(z - 40.0d0))**2
        cz = 1.0d0 - cz + 0.5d0*cz**2
        do 10 jx = 1,mx
          x = (jx - 1.0d0)*dx
          cx = (0.1d0*(x - 10.0d0))**2
          cx = 1.0d0 - cx + 0.5d0*cx**2
          y(1,jx,jz) = 1.0d6*cx*cz
          y(2,jx,jz) = 1.0d12*cx*cz
 10       continue
 20     continue
c
c Write heading, problem parameters, solution parameters.
      write(6,30) mx, mz, mf, rtol, atol
 30   format('Demonstration program for DLSODKR package'//
     1       '2D diurnal kinetics-transport PDE system with 2 species'/
     2       'Spatial mesh is',i3,' by',i3/'Method flag is mf =',i3,
     3       '   Tolerances are rtol =',d8.1,'   atol =',d8.1/
     4       'Left preconditioner uses block-diagonal part of Jacobian'/
     5       'Root function finds stationary points of total ozone,'/
     6       '  i.e. roots of (d/dt)(sum of c2 over all mesh points)'/)
c
c Loop over output points, call DLSODKR, print sample solution values.
      do 70 iout = 1,13
 40     call dlsodkr (fdem, neq, y, t, tout, 1, rtol, atol, 1, istate,
     1      0, rwork, lrw, iwork, liw, jacbd, solbd, mf, gdem, 1, jroot)
        write(6,50) t,iwork(11),iwork(14),rwork(11)
 50     format(/' t =',d10.3,4x,'no. steps =',i5,
     1          '   order =',i2,'   stepsize =',d10.3)
        call c2sum (y, mx, mz, c2tot)
        write(6,60) y(1,1,1), y(1,5,5), y(1,10,10),
     1              y(2,1,1), y(2,5,5), y(2,10,10)
 60     format('   c1 (bot.left/middle/top rt.) =',3d12.3/
     1         '   c2 (bot.left/middle/top rt.) =',3d12.3)
        write(6,62)c2tot,jroot
 62     format('   total c2 =',d15.6,
     1         '   jroot =',i2' (1 = root found, 0 = no root)')
        if (istate .lt. 0) then
          write(6,65)istate
 65       format('DLSODKR returned istate = ',i3)
          go to 75
        endif
        if (istate .eq. 3) then
          istate = 2
          go to 40
        endif
        tout = tout + twohr
 70     continue
c
c Print final statistics.
 75   lenrw = iwork(17)
      leniw = iwork(18)
      nst = iwork(11)
      nsfi = iwork(24)
      nfe = iwork(12)
      nge = iwork(10)
      npe = iwork(13)
      njev = iwork(25)
      nps = iwork(21)
      nni = iwork(19)
      nli = iwork(20)
      avdim = real(nli)/real(nni)
      ncfn = iwork(22)
      ncfl = iwork(23)
      write (6,80) lenrw,leniw,nst,nsfi,nfe,nge,npe,njev,nps,nni,nli,
     1             avdim,ncfn,ncfl
 80   format(//' Final statistics:'/
     1 ' rwork size =',i5,5x,' iwork size =',i4/
     2 ' number of steps        =',i5,5x,'no. fnal. iter. steps  =',i5/
     3 ' number of f evals.     =',i5,5x,'number of g evals.     =',i5/
     4 ' number of prec. evals. =',i5,5x,'number of Jac. evals.  =',i5/
     5 ' number of prec. solves =',i5/
     6 ' number of nonl. iters. =',i5,5x,'number of lin. iters.  =',i5/
     7 ' average Krylov subspace dimension (nli/nni)  =',f8.4/
     8 ' number of conv. failures:  nonlinear =',i3,'  linear =',i3)
      stop
      end

      subroutine fdem (neq, t, y, ydot)
      integer neq,  mx, mz, mm,
     1   iblok, iblok0, idn, iup, ileft, iright, jx, jz
      double precision t, y(2,*), ydot(2,*),
     1   q1, q2, q3, q4, a3, a4, om, c3, dz, hdco, vdco, haco,
     2   c1, c1dn, c1up, c1lt, c1rt, c2, c2dn, c2up, c2lt, c2rt, 
     3   czdn, czup, horad1, hord1, horad2, hord2, qq1, qq2, qq3, qq4, 
     4   rkin1, rkin2, s, vertd1, vertd2, zdn, zup
      common /pcom/ q1,q2,q3,q4,a3,a4,om,c3,dz,hdco,vdco,haco,mx,mz,mm
c
c Set diurnal rate coefficients.
      s = sin(om*t)
      if (s .gt. 0.0d0) then
        q3 = exp(-a3/s)
        q4 = exp(-a4/s)
      else
        q3 = 0.0d0
        q4 = 0.0d0
      endif
c Loop over all grid points.
      do 20 jz = 1,mz
        zdn = 30.0d0 + (jz - 1.5d0)*dz
        zup = zdn + dz
        czdn = vdco*exp(0.2d0*zdn)
        czup = vdco*exp(0.2d0*zup)
        iblok0 = (jz-1)*mx
        idn = -mx
        if (jz .eq. 1) idn = mx
        iup = mx
        if (jz .eq. mz) iup = -mx
        do 10 jx = 1,mx
          iblok = iblok0 + jx
          c1 = y(1,iblok)
          c2 = y(2,iblok)
c Set kinetic rate terms.
          qq1 = q1*c1*c3
          qq2 = q2*c1*c2
          qq3 = q3*c3
          qq4 = q4*c2
          rkin1 = -qq1 - qq2 + 2.0d0*qq3 + qq4
          rkin2 = qq1 - qq2 - qq4
c Set vertical diffusion terms.
          c1dn = y(1,iblok+idn)
          c2dn = y(2,iblok+idn)
          c1up = y(1,iblok+iup)
          c2up = y(2,iblok+iup)
          vertd1 = czup*(c1up - c1) - czdn*(c1 - c1dn)
          vertd2 = czup*(c2up - c2) - czdn*(c2 - c2dn)
c Set horizontal diffusion and advection terms.
          ileft = -1
          if (jx .eq. 1) ileft = 1
          iright = 1
          if (jx .eq. mx) iright = -1
          c1lt = y(1,iblok+ileft)
          c2lt = y(2,iblok+ileft)
          c1rt = y(1,iblok+iright)
          c2rt = y(2,iblok+iright)
          hord1 = hdco*(c1rt - 2.0d0*c1 + c1lt)
          hord2 = hdco*(c2rt - 2.0d0*c2 + c2lt)
          horad1 = haco*(c1rt - c1lt)
          horad2 = haco*(c2rt - c2lt)
c Load all terms into ydot.
          ydot(1,iblok) = vertd1 + hord1 + horad1 + rkin1
          ydot(2,iblok) = vertd2 + hord2 + horad2 + rkin2
 10       continue
 20     continue
      return
      end

      subroutine gdem (neq, t, y, ng, gout)
      integer neq, ng,  mx, mz, mm,
     1   iblok, iblok0, idn, iup, ileft, iright, jx, jz
      double precision t, y(2,*), gout,
     1   q1, q2, q3, q4, a3, a4, om, c3, dz, hdco, vdco, haco,
     2   c1, c2, c2dn, c2up, c2lt, c2rt, c2dot, czdn, czup, horad2,
     3   hord2, qq1, qq2, qq4, rkin2, s, sum, vertd2, zdn, zup
      common /pcom/ q1,q2,q3,q4,a3,a4,om,c3,dz,hdco,vdco,haco,mx,mz,mm
c
c This routine computes the rates for c2 and adds them.
c
c Set diurnal rate coefficient q4.
      s = sin(om*t)
      if (s .gt. 0.0d0) then
        q4 = exp(-a4/s)
      else
        q4 = 0.0d0
      endif
      sum = 0.0d0
c Loop over all grid points.
      do 20 jz = 1,mz
        zdn = 30.0d0 + (jz - 1.5d0)*dz
        zup = zdn + dz
        czdn = vdco*exp(0.2d0*zdn)
        czup = vdco*exp(0.2d0*zup)
        iblok0 = (jz-1)*mx
        idn = -mx
        if (jz .eq. 1) idn = mx
        iup = mx
        if (jz .eq. mz) iup = -mx
        do 10 jx = 1,mx
          iblok = iblok0 + jx
          c1 = y(1,iblok)
          c2 = y(2,iblok)
c Set kinetic rate term for c2.
          qq1 = q1*c1*c3
          qq2 = q2*c1*c2
          qq4 = q4*c2
          rkin2 = qq1 - qq2 - qq4
c Set vertical diffusion terms for c2.
          c2dn = y(2,iblok+idn)
          c2up = y(2,iblok+iup)
          vertd2 = czup*(c2up - c2) - czdn*(c2 - c2dn)
c Set horizontal diffusion and advection terms for c2.
          ileft = -1
          if (jx .eq. 1) ileft = 1
          iright = 1
          if (jx .eq. mx) iright = -1
          c2lt = y(2,iblok+ileft)
          c2rt = y(2,iblok+iright)
          hord2 = hdco*(c2rt - 2.0d0*c2 + c2lt)
          horad2 = haco*(c2rt - c2lt)
c Load all terms into c2dot and sum.
          c2dot = vertd2 + hord2 + horad2 + rkin2
          sum = sum + c2dot
 10       continue
 20     continue
      gout = sum
      return
      end

      subroutine jacbd (f, neq, t, y, ysv, rewt, f0, f1, hl0, jok,
     1    bd, ipbd, ier)
      external f
      integer neq, jok, ipbd(2,*), ier, mx, mz, mm,
     1    iblok, iblok0, jx, jz, lenbd
      double precision t, y(2,*), ysv(neq), rewt(neq), f0(neq), f1(neq),
     1   hl0, bd(2,2,*),
     2   q1, q2, q3, q4, a3, a4, om, c3, dz, hdco, vdco, haco,
     3   c1, c2, czdn, czup, diag, temp, zdn, zup
      common /pcom/ q1,q2,q3,q4,a3,a4,om,c3,dz,hdco,vdco,haco,mx,mz,mm
c
        lenbd = 4*mm
c If jok = 1, copy saved block-diagonal approximate Jacobian into bd.
      if (jok .eq. 1) then
        call dcopy (lenbd, bd(1,1,mm+1), 1, bd, 1)
        go to 30
        endif
c
c If jok = -1, compute and save diagonal Jacobian blocks
c  (using q3 and q4 values computed on last f call).
      do 20 jz = 1,mz
        zdn = 30.0d0 + (jz - 1.5d0)*dz
        zup = zdn + dz
        czdn = vdco*exp(0.2d0*zdn)
        czup = vdco*exp(0.2d0*zup)
        diag = -(czdn + czup + 2.0d0*hdco)
        iblok0 = (jz-1)*mx
        do 10 jx = 1,mx
          iblok = iblok0 + jx
          c1 = y(1,iblok)
          c2 = y(2,iblok)
          bd(1,1,iblok) = (-q1*c3 - q2*c2) + diag
          bd(1,2,iblok) = -q2*c1 + q4
          bd(2,1,iblok) = q1*c3 - q2*c2
          bd(2,2,iblok) = (-q2*c1 - q4) + diag
 10       continue
 20     continue
      call dcopy (lenbd, bd, 1, bd(1,1,mm+1), 1)
c Scale by -hl0, add identity matrix and LU-decompose blocks.
 30   temp = -hl0
      call dscal (lenbd, temp, bd, 1)
      do 40 iblok = 1,mm
        bd(1,1,iblok) = bd(1,1,iblok) + 1.0d0
        bd(2,2,iblok) = bd(2,2,iblok) + 1.0d0
        call dgefa (bd(1,1,iblok), 2, 2, ipbd(1,iblok), ier)
        if (ier .ne. 0) return
 40     continue
      return
      end

      subroutine solbd (neq, t, y, f0, wk, hl0, bd, ipbd, v, lr, ier)
      integer neq, ipbd(2,*), lr, ier,  mx, mz, mm,  i
      double precision t, y(neq), f0(neq), wk(neq), hl0, bd(2,2,*),
     2   v(2,*),  q1, q2, q3, q4, a3, a4, om, c3, dz, hdco, vdco, haco
      common /pcom/ q1,q2,q3,q4,a3,a4,om,c3,dz,hdco,vdco,haco,mx,mz,mm
c Solve the block-diagonal system Px = v using LU factors stored in bd
c and pivot data in ipbd, and return the solution in v.
      ier = 0
      do 10 i = 1,mm
        call dgesl (bd(1,1,i), 2, 2, ipbd(1,i), v(1,i), 0)
 10     continue
      return
      end

      subroutine c2sum (y, mx, mz, c2tot)
      integer mx, mz, jx, jz
      double precision y(2,mx,mz), c2tot, sum
c Sum the c2 values.
      sum = 0.0d0
      do 20 jz = 1,mz
        do 20 jx = 1,mx
 20       sum = sum + y(2,jx,jz)
      c2tot = sum
      return
      end
