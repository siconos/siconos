c-----------------------------------------------------------------------
c Demonstration program for the DLSOIBT package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c This program solves a semi-discretized form of the following system
c Of three PDEs (each similar to a Burgers equation):
c
c   u(i)   =  -(u(1)+u(2)+u(3)) u(i)   +  eta(i) u(i)    (i=1,2,3),
c       t                           x                xx
c
c on the interval  -1 .le. x .le. 1, and with time t .ge. 0.
c The diffusion coefficients are eta(*) = .1, .02, .01.
c The boundary conditions are u(i) = 0 at x = -1 and x = 1 for all i.
c The initial profile for each u(i) is a square wave:
c     u(i) = 0         on 1/2 .lt. abs(x) .le. 1
c     u(i) = amp(i)/2  on abs(x) = 1/2
c     u(i) = amp(i)    on 0 .le. abs(x) .lt. 1/2
c where the amplitudes are amp(*) = .2, .3, .5.
c
c A simplified Galerkin treatment of the spatial variable x is used,
c with piecewise linear basis functions on a uniform mesh of 100
c intervals.  The result is a system of ODEs in the discrete values
c u(i,k) approximating u(i)  (i=1,2,3) at the interior points
c (k = 1,...,99).  The ODEs are:
c
c    .            .        .
c   (u(i,k-1) + 4 u(i,k) + u(i,k+1))/6  =
c
c     -(1/6dx) (c(k-1)dul(i) + 2c(k)(dul(i)+dur(i)) + c(k+1)dur(i))
c
c     + (eta(i)/dx**2) (dur(i) - dul(i))     (i=1,2,3,  k=1,...,99),
c
c where
c     c(j) = u(1,j)+u(2,j)+u(3,j),   dx = .02 = the interval size,
c     dul(i) = u(i,k) - u(i,k-1),   dur(i) = u(i,k+1) - u(i,k).
c Terms involving boundary values (subscripts 0 or 100) are dropped
c from the equations for k = 1 and k = 99 above.
c
c The problem is run for each of the 4 values of mf, and for two values
c of the tolerances.  Output is taken at t = .1, .2, .3, .4.
c Output is on unit lout, set to 6 in a data statement below.
c-----------------------------------------------------------------------
      external res, addabt, jacbt
      integer ncomp, nip, nm1
      integer i, io, istate, itol, iwork, jtol, lout, liw, lrw,
     1   meth, miter, mf, neq, nerr, nint, nout
      double precision eodsq, r6d
      double precision abermx, atol, dx, errfac, eta, hun, one,
     1   rtol, rwork, six, t, tinit, tlast, tout, tols, two, y, ydoti
      dimension eta(3), y(297), ydoti(297), tout(4), tols(2)
      dimension rwork(7447), iwork(317)
c Pass problem parameters in the common block par.
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
c
c Set problem parameters and run parameters
      data eta/0.1d0,0.02d0,0.01d0/, tinit/0.0d0/, tlast/0.4d0/
      data one/1.0d0/, two/2.0d0/, six/6.0d0/, hun/100.0d0/
      data tout/.10d0,.20d0,.30d0,.40d0/
      data lout/6/, nout/4/, lrw/7447/, liw/317/
      data itol/1/, tols/1.0d-3, 1.0d-6/
c
c Set mesh parameters nint, dxc etc.
      nint = 100
      ncomp = 3
      dx = two/nint
      r6d = one/(six*dx)
      do 10 i = 1,ncomp
 10     eodsq(i) = eta(i)/dx**2
      nip = nint - 1
      neq = ncomp*nip
      nm1 = nip - 1
      iwork(1) = ncomp
      iwork(2) = nip
c
      nerr = 0
c
c Set the initial conditions (for output purposes only).
      call setic (nint, ncomp, y)
c
      write (lout,1000)
      write (lout,1100) (eta(i),i=1,ncomp), tinit, tlast, nint,
     1      ncomp, nip, neq
      write (lout,1200)
      call edit (y, ncomp, nip, lout)
c
c The jtol loop is over error tolerances.
      do 200 jtol = 1,2
      rtol = tols(jtol)
      atol = rtol
c
c The meth/miter loops cover 4 values of method flag mf.
      do 100 meth = 1,2
       do 100 miter = 1,2
        mf = 10*meth + miter
c
c Set the initial conditions.
        call setic (nint, ncomp, y)
        t = tinit
        istate = 0
c
        write (lout,1500)  rtol, atol, mf
c
c Loop over output times for each case
        do 80 io = 1,nout
c
          call dlsoibt (res, addabt,jacbt, neq, y, ydoti, t, tout(io),
     1     itol,rtol,atol, 1, istate, 0, rwork,lrw,iwork,liw, mf)
c
          write (lout,2000) t, rwork(11), iwork(14), iwork(11)
          if (io .eq. nout) call edit (y, ncomp, nip, lout)
c
c If istate is not 2 on return, print message and go to next case.
          if (istate .ne. 2) then
            write (lout,4000) mf, t, istate
            nerr = nerr + 1
            go to 100
            endif
c
 80       continue
c
c Print final statistics.
        write (lout,3000) mf, iwork(11), iwork(12), iwork(13),
     1         iwork(17), iwork(18)
c
c Estimate final error and print result.
        call maxerr (y, ncomp, nip, abermx)
        errfac = abermx/tols(jtol)
        if (errfac .lt. hun) then
          write (lout,5000) errfac
        else  
          write (lout,5100) errfac
          nerr = nerr + 1
          endif
 100    continue
 200  continue
c
      write (lout,6000) nerr
      stop
c
 1000 format(/20x,' Demonstration Problem for DLSOIBT'//
     1   10x,'Galerkin method solution of system of 3 PDEs:'//
     2   10x,'  u(i)   =  -(u(1)+u(2)+u(3)) u(i)   +  eta(i) u(i)',
     3   5x,'(i=1,2,3)',/16x,'t',27x,'x',16x,'xx'//
     4   10x,'x interval is -1 to 1,  zero boundary conditions'/
     5   10x,'x discretized using piecewise linear basis functions')
 1100 format(/10x,'Fixed parameters are as follows:'/
     1       13x,'Diffusion coefficients are eta =',3d10.2/
     2       13x,'t0 = ',d12.5/13x,'tlast = ',d12.5/
     3       13x,'Uniform mesh, number of intervals =',i4/
     4       13x,'Block size mb =',i2/13x,'Number of blocks nb =',i4/
     5       13x,'ODE system size neq =',i5//)
c
 1200 format(/'Initial profiles:'/)
c
 1500 format(////90('*')//'Run with rtol =',d9.1,'  atol =',d9.1,
     1       '   mf =',i3///)
c
 2000 format(' At time t =',d12.5,'  current h =',d12.5,
     1       '  current order =',i2,'  current nst =',i5/)
c
 3000 format(//'Final statistics for mf = ',i2,':',
     1       i5,' steps,',i6,' res,',i6,' jacobians,'/
     2       30x,       'rwork size =',i6,',  iwork size =',i6)
c
 4000 format(//20x,'Final time reached for mf = ',i2,
     1       ' was t = ',d12.5/25x,'at which istate = ',i2//)
 5000 format('Final output is correct to within ',d9.2,
     1       '  times local error tolerance. ')
 5100 format('Final output is wrong by ',d9.2,
     1       '  times local error tolerance.')
 6000 format(//90('*')//'Run completed: ',i3,' errors encountered')
c
c end of main program for the DLSOIBT demonstration problem
      end

      subroutine setic (nint, mb, y)
c This routine loads the y array with initial data based on a
c square wave profile for each of the mb PDE variables.
c
      integer nint, mb, i, k, nip, n14, n34
      double precision y,  amp, half, zero
      dimension y(mb,*), amp(3)
      data zero/0.0d0/, half/0.5d0/, amp/0.2d0,0.3d0,0.5d0/
c
      nip = nint - 1
      n14 = nint/4
      n34 = 3*n14
c
      do 15 k = 1,n14-1
        do 10 i = 1,mb
 10       y(i,k) = zero
 15     continue
c
      do 20 i = 1,mb
 20     y(i,n14) = half*amp(i)
c
      do 35 k = n14+1,n34-1
        do 30 i = 1,mb
 30       y(i,k) = amp(i)
 35     continue
c
      do 40 i = 1,mb
 40     y(i,n34) = half*amp(i)
c
      do 55 k = n34+1,nip
        do 50 i = 1,mb
 50       y(i,k) = zero
 55     continue
c
      return
c end of subroutine setic
      end

      subroutine res (n, t, y, v, r, ires)
c This subroutine computes the residual vector
c   r = g(t,y) - A(t,y)*v
c using routines gfun and subav.
c If ires = -1, only g(t,y) is returned in r, since A(t,y) does
c not depend on y.
c No changes need to be made to this routine if nip is changed.
c
      integer ires, n,  ncomp, nip, nm1
      double precision t, y, v, r, r6d, eodsq
      dimension y(n), v(n), r(n)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
c
      call gfun (t, y, r, ncomp)
      if (ires .eq. -1) return
c
      call subav (r, v, ncomp)
c
      return
c end of subroutine res
      end

      subroutine gfun (t, y, g, mb)
c This subroutine computes the right-hand side function g(y,t).
c It uses r6d = 1/(6*dx), eodsq(*) = eta(*)/dx**2, nip,
c and nm1 = nip - 1 from the Common block par.
c
      integer mb,  ncomp, nip, nm1,  i, k
      double precision t, y, g,  r6d, eodsq,  cc, cl, cr, dli, dri, two
      dimension g(mb,*), y(mb,*)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
      data two/2.0d0/
c
c left-most interior point (k = 1)
      cc = y(1,1) + y(2,1) + y(3,1)
      cr = y(1,2) + y(2,2) + y(3,2)
      do 10 i = 1,mb
        dri = y(i,2) - y(i,1)
        g(i,1) = -r6d*(two*cc*y(i,2) + cr*dri)
     1         + eodsq(i)*(dri - y(i,1))
 10     continue
c
c interior points k = 2 to nip-1
      do 20 k = 2,nm1
        cl = y(1,k-1) + y(2,k-1) + y(3,k-1)
        cc = y(1,k) + y(2,k) + y(3,k)
        cr = y(1,k+1) + y(2,k+1) + y(3,k+1)
        do 15 i = 1,mb
          dli = y(i,k) - y(i,k-1)
          dri = y(i,k+1) - y(i,k)
          g(i,k) = -r6d*(cl*dli + two*cc*(dli + dri) + cr*dri)
     1           + eodsq(i)*(dri - dli)
 15       continue
 20     continue
c
c right-most interior point (k = nip)
      cl = y(1,nm1) + y(2,nm1) + y(3,nm1)
      cc = y(1,nip) + y(2,nip) + y(3,nip)
      do 30 i = 1,mb
        dli = y(i,nip) - y(i,nm1)
        g(i,nip) = -r6d*(cl*dli - two*cc*y(i,nm1))
     1           - eodsq(i)*(y(i,nip) + dli)
 30     continue
c
      return
c end of subroutine gfun
      end

      subroutine subav (r, v, mb)
c This routine subtracts the matrix a time the vector v from r,
c in order to form the residual vector, stored in r.
c
      integer mb,  ncomp, nip, nm1,  i, k
      double precision r, v,  r6d, eodsq,  aa1, aa4, four, one, six
      dimension r(mb,*), v(mb,*)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
      data one /1.0d0/, four /4.0d0/, six /6.0d0/
c
      aa1 = one/six
      aa4 = four/six
c
      do 10 i = 1,mb
 10     r(i,1) = r(i,1) - (aa4*v(i,1) + aa1*v(i,2))
c
      do 20 k = 2,nm1
        do 15 i = 1,mb
 15       r(i,k) = r(i,k) - (aa1*v(i,k-1) + aa4*v(i,k) + aa1*v(i,k+1))
 20     continue
c
      do 30 i = 1,mb
 30     r(i,nip) = r(i,nip) - (aa1*v(i,nm1) + aa4*v(i,nip))
c
      return
c end of subroutine subav
      end

      subroutine addabt (n, t, y, mb, nb, pa, pb, pc)
c This subroutine computes the elements of the matrix A,
c and adds them to pa, pb, and pc in the appropriate manner.
c The matrix A is tridiagonal, of order n, with
c nonzero elements (reading across) of 1/6, 4/6, 1/6.
c
      integer n, mb, nb, i, k
      double precision pa, pb, pc, t, y,  aa1, aa4, four, one, six
      dimension y(mb,nb),pa(mb,mb,nb),pb(mb,mb,nb),pc(mb,mb,nb)
      data one/1.0d0/, four/4.0d0/, six/6.0d0/
c
      aa1 = one/six
      aa4 = four/six
      do 50 k = 1,nb
        do 10 i = 1,mb
 10       pa(i,i,k) = pa(i,i,k) + aa4
        if (k .ne. nb) then
          do 20 i = 1,mb
 20         pb(i,i,k) = pb(i,i,k) + aa1
        endif
        if (k .ne. 1) then
          do 30 i = 1,mb
 30         pc(i,i,k) = pc(i,i,k) + aa1
        endif
 50     continue
c
      return
c end of subroutine addabt
      end

      subroutine jacbt (n, t, y, s, mb, nb, pa, pb, pc)
c This subroutine computes the Jacobian dg/dy = d(g-A*s)/dy
c which has block-tridiagonal structure.  The main, upper, and
c lower diagonals are stored in pa, pb, and pc, respectively.
c
      integer n, mb, nb,  ncomp, nip, nm1, i, j, k
      double precision t, y, s, pa, pb, pc,  r6d, eodsq,  cc, cl, cr,
     1   dlj, drj, paij, pbij, pcij, terma, termb, termc, two
      dimension y(mb,nb),s(n),pa(mb,mb,nb),pb(mb,mb,nb),pc(mb,mb,nb)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
      data two/2.0d0/
c
c left-most interior point (k = 1)
      cc = y(1,1) + y(2,1) + y(3,1)
      cr = y(1,2) + y(2,2) + y(3,2)
      terma = r6d*cr
      termb = -r6d*(two*cc + cr)
      do 20 j = 1,mb
        drj = y(j,2) - y(j,1)
        paij = -r6d*two*y(j,2)
        pbij = -r6d*drj
        do 10 i = 1,mb
          pa(i,j,1) = paij
 10       pb(i,j,1) = pbij
        pa(j,j,1) = pa(j,j,1) + terma - two*eodsq(j)
        pb(j,j,1) = pb(j,j,1) + termb + eodsq(j)
 20     continue
c
c interior points k = 2 to nip-1
      do 50 k = 2,nm1
        cl = y(1,k-1) + y(2,k-1) + y(3,k-1)
        cc = y(1,k) + y(2,k) + y(3,k)
        cr = y(1,k+1) + y(2,k+1) + y(3,k+1)
        terma = r6d*(cr - cl)
        termb = -r6d*(two*cc + cr)
        termc = r6d*(two*cc + cl)
        do 40 j = 1,mb
          dlj = y(j,k) - y(j,k-1)
          drj = y(j,k+1) - y(j,k)
          paij = -r6d*two*(dlj + drj)
          pbij = -r6d*drj
          pcij = -r6d*dlj
          do 30 i = 1,mb
            pa(i,j,k) = paij
            pb(i,j,k) = pbij
 30         pc(i,j,k) = pcij
          pa(j,j,k) = pa(j,j,k) + terma - two*eodsq(j)
          pb(j,j,k) = pb(j,j,k) + termb + eodsq(j)
          pc(j,j,k) = pc(j,j,k) + termc + eodsq(j)
 40       continue
 50     continue
c
c right-most interior point (k = nip)
      cl = y(1,nm1) + y(2,nm1) + y(3,nm1)
      cc = y(1,nip) + y(2,nip) + y(3,nip)
      terma = -r6d*cl
      termc = r6d*(two*cc + cl)
      do 70 j = 1,mb
        dlj = y(j,nip) - y(j,nm1)
        paij = r6d*two*y(j,nm1)
        pcij = -r6d*dlj
        do 60 i = 1,mb
          pa(i,j,nip) = paij
 60       pc(i,j,nip) = pcij
        pa(j,j,nip) = pa(j,j,nip) + terma - two*eodsq(j)
        pc(j,j,nip) = pc(j,j,nip) + termc + eodsq(j)
 70     continue
c
      return
c end of subroutine jacbt
      end

      subroutine edit (y, mb, nip, lout)
c This routine prints output.  For each of the mb PDE components, the
c values at the nip points are printed.  All output is on unit lout.
c
      integer mb, nip, lout,  i, k
      double precision y
      dimension y(mb,nip)
c
      do 10 i = 1,mb
 10      write (lout,20) i, (y(i,k),k=1,nip)
c
 20   format(' Values of PDE component i =',i3/15(7d12.4/))
c
      return
c end of subroutine edit
      end

      subroutine maxerr (y, mb, nb, abermx)
c This routine computes the maximum absolute error in the y array,
c as a computed solution at t = 0.4, using data-loaded values for
c accurate answers (from running with much smaller tolerances).
c
      integer mb, nb,  k
      double precision y, abermx,  ae1, ae2, ae3, u1, u2, u3, zero,
     1   u1a, u1b, u1c, u1d, u1e, u1f, u1g,
     2   u2a, u2b, u2c, u2d, u2e, u2f, u2g,
     3   u3a, u3b, u3c, u3d, u3e, u3f, u3g
      dimension y(mb,nb), u1(99), u2(99), u3(99)
      dimension u1a(16),u1b(16),u1c(16),u1d(16),u1e(16),u1f(16),u1g(3),
     2          u2a(16),u2b(16),u2c(16),u2d(16),u2e(16),u2f(16),u2g(3),
     3          u3a(16),u3b(16),u3c(16),u3d(16),u3e(16),u3f(16),u3g(3)
      equivalence (u1a(1),u1(1)), (u1b(1),u1(17)),
     1      (u1c(1),u1(33)), (u1d(1),u1(49)), (u1e(1),u1(65)),
     1      (u1f(1),u1(81)), (u1g(1),u1(97))
      equivalence (u2a(1),u2(1)), (u2b(1),u2(17)),
     1      (u2c(1),u2(33)), (u2d(1),u2(49)), (u2e(1),u2(65)),
     1      (u2f(1),u2(81)), (u2g(1),u2(97))
      equivalence (u3a(1),u3(1)), (u3b(1),u3(17)),
     1      (u3c(1),u3(33)), (u3d(1),u3(49)), (u3e(1),u3(65)),
     1      (u3f(1),u3(81)), (u3g(1),u3(97))
c
      data u1a /
     1  1.70956682d-03, 3.43398445d-03, 5.18783349d-03, 6.98515842d-03,
     1  8.83921016d-03, 1.07622016d-02, 1.27650806d-02, 1.48573251d-02,
     1  1.70467655d-02, 1.93394396d-02, 2.17394852d-02, 2.42490773d-02,
     1  2.68684152d-02, 2.95957660d-02, 3.24275691d-02, 3.53586054d-02/
      data u1b /
     1  3.83822285d-02, 4.14906520d-02, 4.46752791d-02, 4.79270545d-02,
     1  5.12368132d-02, 5.45956048d-02, 5.79949684d-02, 6.14271460d-02,
     1  6.48852271d-02, 6.83632267d-02, 7.18561029d-02, 7.53597274d-02,
     1  7.88708192d-02, 8.23868545d-02, 8.59059616d-02, 8.94268082d-02/
      data u1c /
     1  9.29484864d-02, 9.64703968d-02, 9.99921344d-02, 1.03513375d-01,
     1  1.07033760d-01, 1.10552783d-01, 1.14069668d-01, 1.17583246d-01,
     1  1.21091827d-01, 1.24593066d-01, 1.28083828d-01, 1.31560049d-01,
     1  1.35016617d-01, 1.38447256d-01, 1.41844451d-01, 1.45199401d-01/
      data u1d /
     1  1.48502033d-01, 1.51741065d-01, 1.54904135d-01, 1.57977973d-01,
     1  1.60948623d-01, 1.63801670d-01, 1.66522463d-01, 1.69096305d-01,
     1  1.71508595d-01, 1.73744902d-01, 1.75790974d-01, 1.77632682d-01,
     1  1.79255895d-01, 1.80646319d-01, 1.81789276d-01, 1.82669470d-01/
      data u1e /
     1  1.83270725d-01, 1.83575716d-01, 1.83565712d-01, 1.83220322d-01,
     1  1.82517279d-01, 1.81432251d-01, 1.79938706d-01, 1.78007835d-01,
     1  1.75608540d-01, 1.72707519d-01, 1.69269456d-01, 1.65257378d-01,
     1  1.60633244d-01, 1.55358941d-01, 1.49398029d-01, 1.42718981d-01/
      data u1f /
     1  1.35301474d-01, 1.27148627d-01, 1.18308730d-01, 1.08905085d-01,
     1  9.91559295d-02, 8.93515884d-02, 7.97824293d-02, 7.06663514d-02,
     1  6.21244732d-02, 5.41994827d-02, 4.68848207d-02, 4.01465202d-02,
     1  3.39357642d-02, 2.81954415d-02, 2.28635569d-02, 1.78750916d-02/
      data u1g / 1.31630892d-02, 8.65933391d-03, 4.29480447d-03/
      data u2a /
     1  7.17416019d-06, 1.70782645d-05, 3.31245126d-05, 6.01588363d-05,
     1  1.05339286d-04, 1.79174771d-04, 2.96719122d-04, 4.78862606d-04,
     1  7.53598916d-04, 1.15707860d-03, 1.73420412d-03, 2.53849668d-03,
     1  3.63099110d-03, 5.07800919d-03, 6.94782549d-03, 9.30645443d-03/
      data u2b /
     1  1.22130079d-02, 1.57152366d-02, 1.98459102d-02, 2.46205841d-02,
     1  3.00370492d-02, 3.60764461d-02, 4.27057301d-02, 4.98809820d-02,
     1  5.75510102d-02, 6.56607602d-02, 7.41541974d-02, 8.29764928d-02,
     1  9.20754824d-02, 1.01402468d-01, 1.10912474d-01, 1.20564094d-01/
      data u2c /
     1  1.30319039d-01, 1.40141489d-01, 1.49997326d-01, 1.59853293d-01,
     1  1.69676126d-01, 1.79431680d-01, 1.89084097d-01, 1.98595037d-01,
     1  2.07923034d-01, 2.17023055d-01, 2.25846345d-01, 2.34340694d-01,
     1  2.42451240d-01, 2.50121934d-01, 2.57297724d-01, 2.63927433d-01/
      data u2d /
     1  2.69967170d-01, 2.75383917d-01, 2.80158840d-01, 2.84289739d-01,
     1  2.87792167d-01, 2.90698875d-01, 2.93057586d-01, 2.94927384d-01,
     1  2.96374262d-01, 2.97466488d-01, 2.98270390d-01, 2.98847025d-01,
     1  2.99249945d-01, 2.99524080d-01, 2.99705593d-01, 2.99822450d-01/
      data u2e /
     1  2.99895431d-01, 2.99939301d-01, 2.99963931d-01, 2.99975129d-01,
     1  2.99974996d-01, 2.99961526d-01, 2.99927041d-01, 2.99854809d-01,
     1  2.99712769d-01, 2.99442742d-01, 2.98942676d-01, 2.98038511d-01,
     1  2.96441259d-01, 2.93684573d-01, 2.89040478d-01, 2.81421884d-01/
      data u2f /
     1  2.69315148d-01, 2.50874185d-01, 2.24457680d-01, 1.89885662d-01,
     1  1.49894358d-01, 1.09927672d-01, 7.54041273d-02, 4.90259517d-02,
     1  3.06080023d-02, 1.85165524d-02, 1.09104125d-02, 6.27726960d-03,
     1  3.53002680d-03, 1.94049735d-03, 1.04218859d-03, 5.45964314d-04/
      data u2g / 2.77379128d-04, 1.33343739d-04, 5.32660444d-05/
      data u3a /
     1  1.86765383d-10, 1.96772458d-09, 1.19111389d-08, 5.54964761d-08,
     1  2.18340713d-07, 7.55899524d-07, 2.35604385d-06, 6.70801745d-06,
     1  1.76224112d-05, 4.30351929d-05, 9.82592148d-05, 2.10736217d-04,
     1  4.26209304d-04, 8.15657041d-04, 1.48160943d-03, 2.56186555d-03/
      data u3b /
     1  4.22851247d-03, 6.68078970d-03, 1.01317466d-02, 1.47903961d-02,
     1  2.08424987d-02, 2.84336008d-02, 3.76573037d-02, 4.85502549d-02,
     1  6.10936693d-02, 7.52198901d-02, 9.08218891d-02, 1.07763660d-01,
     1  1.25889931d-01, 1.45034247d-01, 1.65025016d-01, 1.85689556d-01/
      data u3c /
     1  2.06856371d-01, 2.28356037d-01, 2.50021072d-01, 2.71685149d-01,
     1  2.93181998d-01, 3.14344301d-01, 3.35002907d-01, 3.54986687d-01,
     1  3.74123404d-01, 3.92241969d-01, 4.09176451d-01, 4.24772089d-01,
     1  4.38893320d-01, 4.51433444d-01, 4.62324969d-01, 4.71549073d-01/
      data u3d /
     1  4.79142163d-01, 4.85197409d-01, 4.89859810d-01, 4.93314543d-01,
     1  4.95770115d-01, 4.97439231d-01, 4.98520996d-01, 4.99187563d-01,
     1  4.99576941d-01, 4.99791928d-01, 4.99903753d-01, 4.99958343d-01,
     1  4.99983239d-01, 4.99993785d-01, 4.99997902d-01, 4.99999367d-01/
      data u3e /
     1  4.99999835d-01, 4.99999965d-01, 4.99999995d-01, 5.00000000d-01,
     1  5.00000000d-01, 4.99999997d-01, 4.99999976d-01, 4.99999863d-01,
     1  4.99999315d-01, 4.99996914d-01, 4.99987300d-01, 4.99951740d-01,
     1  4.99829328d-01, 4.99435130d-01, 4.98245007d-01, 4.94883400d-01/
      data u3f /
     1  4.86081966d-01, 4.65174923d-01, 4.21856650d-01, 3.47885738d-01,
     1  2.49649938d-01, 1.51648615d-01, 7.80173239d-02, 3.47983164d-02,
     1  1.38686441d-02, 5.05765688d-03, 1.71052539d-03, 5.38966324d-04,
     1  1.57923694d-04, 4.27352191d-05, 1.05512005d-05, 2.33068621d-06/
      data u3g / 4.45404604d-07, 6.88336884d-08, 7.23875975d-09/
      data zero/0.0d0/
c
      abermx = zero
      do 10 k = 1,99
        ae1 = abs(y(1,k) - u1(k))
        ae2 = abs(y(2,k) - u2(k))
        ae3 = abs(y(3,k) - u3(k))
        abermx = max(abermx, ae1, ae2, ae3)
 10     continue
c
      return
c end of subroutine maxerr
      end
