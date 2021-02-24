c-----------------------------------------------------------------------
c Demonstration program for the DLSODAR package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c The DLSODAR package is used to solve two simple problems,
c one nonstiff and one intermittently stiff.
c If the errors are too large, or other difficulty occurs,
c a warning message is printed.  All output is on unit lout = 6.
c-----------------------------------------------------------------------
      external f1, gr1, f2, jac2, gr2
      integer iopt, iout, istate, itask, itol, iwork, jroot, jt,
     1   kroot, leniw, lenrw, liw, lrw, lout, neq, nerr, ng,
     2   nfe, nfea, nge, nje, nst
      double precision atol, er, ero, errt, rtol, rwork,
     1   t, tout, tzero, y, yt
      dimension y(2), atol(2), rwork(57), iwork(22), jroot(2)
      data lout/6/
c
      nerr = 0
c-----------------------------------------------------------------------
c First problem.
c The initial value problem is:
c   dy/dt = ((2*log(y) + 8)/t - 5)*y,  y(1) = 1,  1 .le. t .le. 6
c The solution is  y(t) = exp(-t**2 + 5*t - 4)
c The two root functions are:
c   g1 = ((2*log(y)+8)/t - 5)*y (= dy/dt)  (with root at t = 2.5),
c   g2 = log(y) - 2.2491  (with roots at t = 2.47 and 2.53)
c-----------------------------------------------------------------------
c Set all input parameters and print heading.
      neq = 1
      y(1) = 1.0d0
      t = 1.0d0
      tout = 2.0d0
      itol = 1
      rtol = 1.0d-6
      atol(1) = 1.0d-6
      itask = 1
      istate = 1
      iopt = 0
      lrw = 44
      liw = 21
      jt = 2
      ng = 2
      write (lout,110) itol,rtol,atol(1),jt
 110  format(/' Demonstration program for DLSODAR package'////
     1  ' First problem'///
     2  ' Problem is  dy/dt = ((2*log(y)+8)/t - 5)*y,  y(1) = 1'//
     3  ' Solution is  y(t) = exp(-t**2 + 5*t - 4)'//
     4  ' Root functions are:'/
     5  10x,' g1 = dy/dt  (root at t = 2.5)'/
     6  10x,' g2 = log(y) - 2.2491  (roots at t = 2.47 and t = 2.53)'//
     7  ' itol =',i3,'   rtol =',d10.1,'   atol =',d10.1//
     8  ' jt =',i3///)
c
c Call DLSODAR in loop over tout values 2,3,4,5,6.
      ero = 0.0d0
      do 180 iout = 1,5
 120    continue
        call dlsodar(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jdum,jt,gr1,ng,jroot)
c
c Print y and error in y, and print warning if error too large.
        yt = exp(-t*t + 5.0d0*t - 4.0d0)
        er = y(1) - yt
        write (lout,130) t,y(1),er
 130    format(' At t =',d15.7,5x,'y =',d15.7,5x,'error =',d12.4)
        if (istate .lt. 0) go to 185
        er = abs(er)/(rtol*abs(y(1)) + atol(1))
        ero = max(ero,er)
        if (er .gt. 1000.0d0) then
          write (lout,140)
 140      format(//' Warning: error exceeds 1000 * tolerance'//)
          nerr = nerr + 1
          endif
        if (istate .ne. 3) go to 175
c
c If a root was found, write results and check root location.
c Then reset istate to 2 and return to DLSODAR call.
        write (lout,150) t,jroot(1),jroot(2)
 150    format(/' Root found at t =',d15.7,5x,'jroot =',2i5)
        if (jroot(1) .eq. 1) errt = t - 2.5d0
        if (jroot(2) .eq. 1 .and. t .le. 2.5d0) errt = t - 2.47d0
        if (jroot(2) .eq. 1 .and. t .gt. 2.5d0) errt = t - 2.53d0
        write (lout,160) errt
 160    format(' Error in t location of root is',d12.4/)
        if (abs(errt) .gt. 1.0d-3) then
          write (lout,170)
 170      format(//' Warning: root error exceeds 1.0d-3'//)
          nerr = nerr + 1
          endif
        istate = 2
        go to 120
c
c If no root found, increment tout and loop back.
 175    tout = tout + 1.0d0
 180    continue
c
c Problem complete.  Print final statistics.
 185  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      nge = iwork(10)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 2) nfea = nfe - neq*nje
      write (lout,190) lenrw,leniw,nst,nfe,nfea,nje,nge,ero
 190  format(//' Final statistics for this run:'/
     1  ' rwork size =',i4,'   iwork size =',i4/
     2  ' number of steps =',i5/
     3  ' number of f-s   =',i5/
     4  ' (excluding j-s) =',i5/
     5  ' number of j-s   =',i5/
     6  ' number of g-s   =',i5/
     7  ' error overrun =',d10.2)
c
c-----------------------------------------------------------------------
c Second problem (Van der Pol oscillator).
c The initial value problem is (after reduction of 2nd order ODE):
c   dy1/dt = y2,  dy2/dt = 100*(1 - y1**2)*y2 - y1,
c   y1(0) = 2,  y2(0) = 0,  0 .le. t .le. 200
c The root function is  g = y1.
c An analytic solution is not known, but the zeros of y1 are known
c to 15 figures for purposes of checking the accuracy.
c-----------------------------------------------------------------------
c Set tolerance parameters and print heading.
      itol = 2
      rtol = 1.0d-6
      atol(1) = 1.0d-6
      atol(2) = 1.0d-4
      write (lout,200) itol,rtol,atol(1),atol(2)
 200  format(////80('*')//' Second problem (Van der Pol oscillator)'//
     1  ' Problem is dy1/dt = y2,  dy2/dt = 100*(1-y1**2)*y2 - y1'/
     2  '            y1(0) = 2,  y2(0) = 0'//
     3  ' Root function is  g = y1'//
     4  ' itol =',i3,'   rtol =',d10.1,'   atol =',2d10.1)
c
c Loop over jt = 1, 2.  Set remaining parameters and print jt.
      do 290 jt = 1,2
      neq = 2
      y(1) = 2.0d0
      y(2) = 0.0d0
      t = 0.0d0
      tout = 20.0d0
      itask = 1
      istate = 1
      iopt = 0
      lrw = 57
      liw = 22
      ng = 1
      write (lout,210) jt
 210  format(///' Solution with jt =',i2//)
c
c Call DLSODAR in loop over tout values 20,40,...,200.
      do 270 iout = 1,10
 220    continue
        call dlsodar(f2,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac2,jt,gr2,ng,jroot)
c
c Print y1 and y2.
        write (lout,230) t,y(1),y(2)
 230    format(' At t =',d15.7,5x,'y1 =',d15.7,5x,'y2 =',d15.7)
        if (istate .lt. 0) go to 275
        if (istate .ne. 3) go to 265
c
c If a root was found, write results and check root location.
c Then reset istate to 2 and return to DLSODAR call.
        write (lout,240) t
 240    format(/' Root found at t =',d15.7)
        kroot = int(t/81.2d0 + 0.5d0)
        tzero = 81.17237787055d0 + (kroot-1)*81.41853556212d0
        errt = t - tzero
        write (lout,250) errt
 250    format(' Error in t location of root is',d12.4//)
        if (abs(errt) .gt. 1.0d-1) then
          write (lout,260)
 260      format(//' Warning: root error exceeds 1.0d-1'//)
          nerr = nerr + 1
          endif
        istate = 2
        go to 220
c
c If no root found, increment tout and loop back.
 265    tout = tout + 20.0d0
 270    continue
c
c Problem complete.  Print final statistics.
 275  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      nge = iwork(10)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 2) nfea = nfe - neq*nje
      write (lout,280) lenrw,leniw,nst,nfe,nfea,nje,nge
 280  format(//' Final statistics for this run:'/
     1  '  rwork size =',i4,'   iwork size =',i4/
     2  '  number of steps =',i5/
     3  '  number of f-s   =',i5/
     4  '  (excluding j-s) =',i5/
     5  '  number of j-s   =',i5/
     6  '  number of g-s   =',i5)
 290  continue
c
      write (lout,300) nerr
 300  format(///' Total number of errors encountered =',i3)
      stop
      end

      subroutine f1 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(1), ydot(1)
      ydot(1) = ((2.0d0*log(y(1)) + 8.0d0)/t - 5.0d0)*y(1)
      return
      end

      subroutine gr1 (neq, t, y, ng, groot)
      integer neq, ng
      double precision t, y, groot
      dimension y(1), groot(2)
      groot(1) = ((2.0d0*log(y(1)) + 8.0d0)/t - 5.0d0)*y(1)
      groot(2) = log(y(1)) - 2.2491d0
      return
      end

      subroutine f2 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(2), ydot(2)
      ydot(1) = y(2)
      ydot(2) = 100.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
      return
      end

      subroutine jac2 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y, pd
      dimension y(2), pd(nrowpd,2)
      pd(1,1) = 0.0d0
      pd(1,2) = 1.0d0
      pd(2,1) = -200.0d0*y(1)*y(2) - 1.0d0
      pd(2,2) = 100.0d0*(1.0d0 - y(1)*y(1))
      return
      end

      subroutine gr2 (neq, t, y, ng, groot)
      integer neq, ng
      double precision t, y, groot
      dimension y(2), groot(1)
      groot(1) = y(1)
      return
      end
