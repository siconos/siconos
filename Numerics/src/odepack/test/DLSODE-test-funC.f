C Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2006.
C Siconos is a program dedicated to modeling, simulation and control
C of non smooth dynamical systems.	
C Siconos is a free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C Siconos is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with Siconos; if not, write to the Free Software
C Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
C
C Contact: Vincent ACARY vincent.acary@inrialpes.fr 
C	
c-----------------------------------------------------------------------
c Demonstration program for the DLSODE package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c The package is used to solve two simple problems,
c one with a full Jacobian, the other with a banded Jacobian,
c with all 8 of the appropriate values of mf in each case.
c If the errors are too large, or other difficulty occurs,
c a warning message is printed.  All output is on unit lout = 6.
c-----------------------------------------------------------------------
      external f1, jac1
      integer i, iopar, iopt, iout, istate, itask, itol, iwork,
     1   leniw, lenrw, liw, lout, lrw, mband, meth, mf, miter,
     2   ml, mu, neq, nerr, nfe, nfea, nje, nout, nqu, nst
      double precision atol, dtout, er, erm, ero, hu, rtol, rwork, t,
     1   tout, tout1, y
      dimension y(25), rwork(697), iwork(45)
      data lout/6/, tout1/1.39283880203d0/, dtout/2.214773875d0/
c
      nerr = 0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-6
      lrw = 697
      liw = 45
      iopt = 0
c
c First problem
c
      neq = 2
      nout = 4
      write (lout,110) neq,itol,rtol,atol
 110  format(/' Demonstration program for DLSODE package'///
     1  ' Problem 1:  Van der Pol oscillator:'/
     2  '  xdotdot - 3*(1 - x**2)*xdot + x = 0, ',
     3  '   x(0) = 2, xdot(0) = 0'/
     4  ' neq =',i2/
     5  ' itol =',i3,'   rtol =',d10.1,'   atol =',d10.1//)
c
      do 195 meth = 1,2
      do 190 miter = 0,3
      mf = 10*meth + miter
      write (lout,120) mf
 120  format(///' Solution with mf =',i3//
     1     5x,'t               x               xdot       nq      h'//)
      t = 0.0d0
      y(1) = 2.0d0
      y(2) = 0.0d0
      itask = 1
      istate = 1
      tout = tout1
      ero = 0.0d0
      do 170 iout = 1,nout
        call dlsode(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac1,mf)
        hu = rwork(11)
        nqu = iwork(14)
        write (lout,140) t,y(1),y(2),nqu,hu
 140    format(d15.5,d16.5,d14.3,i5,d14.3)
        if (istate .lt. 0) go to 175
        iopar = iout - 2*(iout/2)
        if (iopar .ne. 0) go to 170
        er = abs(y(1))/atol
        ero = max(ero,er)
        if (er .gt. 1000.0d0) then
          write (lout,150)
 150      format(//' Warning: error exceeds 1000 * tolerance'//)
          nerr = nerr + 1
        endif
 170    tout = tout + dtout
 175  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (miter .eq. 2) nfea = nfe - neq*nje
      if (miter .eq. 3) nfea = nfe - nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 180  format(//' Final statistics for this run:'/
     1  ' rwork size =',i4,'   iwork size =',i4/
     2  ' number of steps =',i5/
     3  ' number of f-s   =',i5/
     4  ' (excluding J-s) =',i5/
     5  ' number of J-s   =',i5/
     6  ' error overrun =',d10.2)
 190  continue
 195  continue

      write (lout,300) nerr
 300  format(////' Number of errors encountered =',i3)
      stop
      end

c$$$      subroutine f1 (neq, t, y, ydot)
c$$$      integer neq
c$$$      double precision t, y, ydot
c$$$      dimension y(neq), ydot(neq)
c$$$      ydot(1) = y(2)
c$$$      ydot(2) = 3.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
c$$$      return
c$$$      end
c$$$
c$$$      subroutine jac1 (neq, t, y, ml, mu, pd, nrowpd)
c$$$      integer neq, ml, mu, nrowpd
c$$$      double precision t, y, pd
c$$$      dimension y(neq), pd(nrowpd,neq)
c$$$      pd(1,1) = 0.0d0
c$$$      pd(1,2) = 1.0d0
c$$$      pd(2,1) = -6.0d0*y(1)*y(2) - 1.0d0
c$$$      pd(2,2) = 3.0d0*(1.0d0 - y(1)*y(1))
c$$$      return
c$$$      end
