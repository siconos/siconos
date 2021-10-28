c-----------------------------------------------------------------------
c Demonstration program for the DLSODA package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c The package is used to solve two simple problems,
c one with a full Jacobian, the other with a banded Jacobian,
c with the 2 appropriate values of jt in each case.
c If the errors are too large, or other difficulty occurs,
c a warning message is printed.  All output is on unit lout = 6.
c-----------------------------------------------------------------------
      external f1, jac1, f2, jac2
      integer i, iopar, iopt, iout, istate, itask, itol, iwork,
     1   jt, leniw, lenrw, liw, lout, lrw, mband, mused,
     2   ml, mu, neq, nerr, nfe, nfea, nje, nout, nqu, nst
      double precision atol, dtout, dtout0, dtout1, er, erm, ero, hu,
     1     rtol, rwork, t, tout, tout1, tsw, y
      dimension y(25), rwork(522), iwork(45)
      data lout/6/, tout1/16.921743d0/, dtout/17.341162d0/
c
      nerr = 0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-8
      lrw = 522
      liw = 45
      iopt = 0
c
c First problem
c
      neq = 2
      nout = 4
      write (lout,110) neq,itol,rtol,atol
 110  format(/'Demonstration program for DLSODA package'////
     1  ' Problem 1:   Van der Pol oscillator:'/
     2  '              xdotdot - 20*(1 - x**2)*xdot + x = 0, ',
     3  '   x(0) = 2, xdot(0) = 0'/' neq =',i2/
     4  ' itol =',i3,'   rtol =',d10.1,'   atol =',d10.1//)
c
      do 190 jt = 1,2
      write (lout,120) jt
 120  format(//' Solution with jt =',i3//
     1       '  t               x               xdot       meth',
     2       '   nq     h           tsw'//)
      t = 0.0d0
      y(1) = 2.0d0
      y(2) = 0.0d0
      itask = 1
      istate = 1
      dtout0 = 0.5d0*tout1
      dtout1 = 0.5d0*dtout
      tout = dtout0
      ero = 0.0d0
      do 170 iout = 1,nout
        call dlsoda(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1              iopt,rwork,lrw,iwork,liw,jac1,jt)
        hu = rwork(11)
        tsw = rwork(15)
        nqu = iwork(14)
        mused = iwork(19)
        write (lout,140) t,y(1),y(2),mused,nqu,hu,tsw
 140    format(d12.5,d16.5,d14.3,2i6,2d13.3)
        if (istate .lt. 0) go to 175
        iopar = iout - 2*(iout/2)
        if (iopar .ne. 0) go to 160
        er = abs(y(1))
        ero = max(ero,er)
        if (er .gt. 1.0d-2) then
          write (lout,150)
 150      format(//' Warning: value at root exceeds 1.0d-2'//)
          nerr = nerr + 1
        endif
 160    if (iout .eq. 1) tout = tout + dtout0
        if (iout .gt. 1) tout = tout + dtout1
 170    continue
 175  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 2) nfea = nfe - neq*nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 180  format(//' Final statistics for this run:'/
     1  ' rwork size =',i4,'   iwork size =',i4/
     2  ' number of steps =',i5/
     3  ' number of f-s   =',i5/
     4  ' (excluding J-s) =',i5/
     5  ' number of J-s   =',i5/
     6  ' max. error at root =',d10.2)
 190  continue
c
c Second problem
c
      neq = 25
      ml = 5
      mu = 0
      iwork(1) = ml
      iwork(2) = mu
      mband = ml + mu + 1
      atol = 1.0d-6
      nout = 5
      write (lout,210) neq,ml,mu,itol,rtol,atol
 210  format(///80('-')///
     1  ' Problem 2: ydot = A * y , where',
     2  '  A is a banded lower triangular matrix'/
     2  '            derived from 2-D advection PDE'/
     3  ' neq =',i3,'   ml =',i2,'   mu =',i2/
     4  ' itol =',i3,'   rtol =',d10.1,'   atol =',d10.1//)
      do 290 jt = 4,5
      write (lout,220) jt
 220  format(//' Solution with jt =',i3//
     1       '     t             max.err.     meth   ',
     2       'nq      h            tsw'//)
      t = 0.0d0
      do 230 i = 2,neq
 230    y(i) = 0.0d0
      y(1) = 1.0d0
      itask = 1
      istate = 1
      tout = 0.01d0
      ero = 0.0d0
      do 270 iout = 1,nout
        call dlsoda(f2,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1              iopt,rwork,lrw,iwork,liw,jac2,jt)
        call edit2(y,t,erm)
        hu = rwork(11)
        tsw = rwork(15)
        nqu = iwork(14)
        mused = iwork(19)
        write (lout,240) t,erm,mused,nqu,hu,tsw
 240    format(d15.5,d14.3,2i6,2d14.3)
        if (istate .lt. 0) go to 275
        er = erm/atol
        ero = max(ero,er)
        if (er .gt. 1000.0d0) then
          write (lout,150)
          nerr = nerr + 1
        endif
 270    tout = tout*10.0d0
 275  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 5) nfea = nfe - mband*nje
      write (lout,280) lenrw,leniw,nst,nfe,nfea,nje,ero
 280  format(//' Final statistics for this run:'/
     1  ' rwork size =',i4,'   iwork size =',i4/
     2  ' number of steps =',i5/
     3  ' number of f-s   =',i5/
     4  ' (excluding J-s) =',i5/
     5  ' number of J-s   =',i5/
     6  ' error overrun =',d10.2)
 290  continue
      write (lout,300) nerr
 300  format(///' Number of errors encountered =',i3)
      stop
      end

      subroutine f1 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(neq), ydot(neq)
      ydot(1) = y(2)
      ydot(2) = 20.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
      return
      end

      subroutine jac1 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y, pd
      dimension y(neq), pd(nrowpd,neq)
      pd(1,1) = 0.0d0
      pd(1,2) = 1.0d0
      pd(2,1) = -40.0d0*y(1)*y(2) - 1.0d0
      pd(2,2) = 20.0d0*(1.0d0 - y(1)*y(1))
      return
      end

      subroutine f2 (neq, t, y, ydot)
      integer neq, i, j, k, ng
      double precision t, y, ydot, alph1, alph2, d
      dimension y(neq), ydot(neq)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      do 10 j = 1,ng
      do 10 i = 1,ng
        k = i + (j - 1)*ng
        d = -2.0d0*y(k)
        if (i .ne. 1) d = d + y(k-1)*alph1
        if (j .ne. 1) d = d + y(k-ng)*alph2
 10     ydot(k) = d
      return
      end

      subroutine jac2 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd, j, mband, mu1, mu2, ng
      double precision t, y, pd, alph1, alph2
      dimension y(neq), pd(nrowpd,neq)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      mband = ml + mu + 1
      mu1 = mu + 1
      mu2 = mu + 2
      do 10 j = 1,neq
        pd(mu1,j) = -2.0d0
        pd(mu2,j) = alph1
 10     pd(mband,j) = alph2
      do 20 j = ng,neq,ng
 20     pd(mu2,j) = 0.0d0
      return
      end

      subroutine edit2 (y, t, erm)
      integer i, j, k, ng
      double precision y, t, erm, alph1, alph2, a1, a2, er, ex, yt
      dimension y(*)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      erm = 0.0d0
      if (t .eq. 0.0d0) return
      ex = 0.0d0
      if (t .le. 30.0d0) ex = exp(-2.0d0*t)
      a2 = 1.0d0
      do 60 j = 1,ng
        a1 = 1.0d0
        do 50 i = 1,ng
          k = i + (j - 1)*ng
          yt = t**(i+j-2)*ex*a1*a2
          er = abs(y(k)-yt)
          erm = max(erm,er)
          a1 = a1*alph1/i
 50       continue
        a2 = a2*alph2/j
 60     continue
      return
      end
