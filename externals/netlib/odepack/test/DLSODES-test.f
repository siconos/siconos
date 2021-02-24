c-----------------------------------------------------------------------
c Demonstration program for the DLSODES package.
c This is the version of 14 June 2001.
c
c This version is in double precision.
c
c The package is used for each of the relevant values of mf to solve
c the problem ydot = A * y, where A is the 9 by 9 sparse matrix
c
c               -4  1     1
c                1 -4  1     1
c                   1 -4        1
c                        -4  1     1
c       A =               1 -4  1     1
c                            1 -4        1
c                                 -4  1
c                                  1 -4  1
c                                     1 -4
c
c The initial conditions are  y(0) = (1, 2, 3, ..., 9).
c Output is printed at t = 1, 2, and 3.
c Each case is solved first with nominal (large) values of lrw and liw,
c and then with values given by lenrw and leniw (optional outputs)
c on the first run, as a check on these computed work array lengths.
c If the errors are too large, or other difficulty occurs,
c a warning message is printed.
c All output is on unit lout, which is data-loaded to 6 below.
c-----------------------------------------------------------------------
      external fdem, jdem
      integer i, ia, igrid, iopt, iout, irun, istate, itask, itol,
     1  iwork, j, ja, k, l, leniw, lenrw, liw, lout, lrw,
     2  m, meth, mf, miter, moss, neq, nerr, nfe, nfea,
     3  ngp, nje, nlu, nnz, nout, nqu, nst, nzl, nzu
      double precision atol, erm, ero, hu, rtol, rwork, t, tout, y
      dimension y(9), ia(10), ja(50), iwork(90), rwork(1000)
      equivalence (ia(1),iwork(31)), (ja(1),iwork(41))
      data lout/6/
c
c Write heading and set fixed parameters.
      write(lout,10)
 10   format(/'Demonstration problem for the DLSODES package'//)
      nerr = 0
      igrid = 3
      neq = igrid**2
      t = 0.0d0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-5
      itask = 1
      iopt = 0
      do 20 i = 1,neq
 20     y(i) = i
      ia(1) = 1
      k = 1
      do 60 m = 1,igrid
        do 50 l = 1,igrid
          j = l + (m - 1)*igrid
          if (m .gt. 1) then
            ja(k) = j - igrid
            k = k + 1
          endif
 30       if (l .gt. 1) then
            ja(k) = j - 1
            k = k + 1
          endif
 35       ja(k) = j
          k = k + 1
          if (l .lt. igrid) then
            ja(k) = j + 1
            k = k + 1
          endif
 40       ia(j+1) = k
 50       continue
 60     continue
      write (lout,80)neq,t,rtol,atol,(y(i),i=1,neq)
 80   format(' neq =',i4,5x,'t0 =',f4.1,5x,'rtol =',d12.3,5x,
     1   'atol =',d12.3//' Initial y vector =  ',9f5.1)
c
c Loop over all relevant values of mf.
      do 193 moss = 0,2
      do 192 meth = 1,2
      do 191 miter = 0,3
      if ( (miter.eq.0 .or. miter.eq.3) .and. moss.ne.0) go to 191
      mf = 100*moss + 10*meth + miter
      write (lout,100)
 100  format(//80('*'))
c First run: nominal work array lengths, 3 output points.
      irun = 1
      lrw = 1000
      liw = 90
      nout = 3
 110  continue
      write (lout,120)mf,lrw,liw
 120  format(//'Run with mf =',i4,'.',5x,
     1       'Input work lengths lrw, liw =',2i6/)
      do 125 i = 1,neq
 125    y(i) = i
      t = 0.0d0
      tout = 1.0d0
      istate = 1
      ero = 0.0d0
c Loop over output points.  Do output and accuracy check at each.
      do 170 iout = 1,nout
        call dlsodes (fdem, neq, y, t, tout, itol, rtol, atol, itask,
     1                istate, iopt, rwork, lrw, iwork, liw, jdem, mf)
        nst = iwork(11)
        hu = rwork(11)
        nqu = iwork(14)
        call edit (y, iout, erm)
        write(lout,140)t,nst,hu,nqu,erm,(y(i),i=1,neq)
 140    format('At t =',f5.1,3x,'nst =',i4,3x,'hu =',d12.3,3x,
     1    'nqu =',i3,3x,' max. err. =',d11.3/
     2    '  y array =    ',4d15.6/5d15.6)
        if (istate .lt. 0) go to 175
        erm = erm/atol
        ero = max(ero,erm)
        if (erm .gt. 100.0d0) then
          write (lout,150)
 150      format(//' Warning: error exceeds 100 * tolerance'//)
          nerr = nerr + 1
        endif
        tout = tout + 1.0d0
 170    continue
 175  continue
      if (istate .lt. 0) nerr = nerr + 1
      if (irun .eq. 2) go to 191
c Print final statistics (first run only)
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nnz = iwork(19)
      ngp = iwork(20)
      nlu = iwork(21)
      nzl = iwork(25)
      nzu = iwork(26)
      nfea = nfe
      if (miter .eq. 2) nfea = nfe - ngp*nje
      if (miter .eq. 3) nfea = nfe - nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 180  format(/'Final statistics for this run:'/
     1  ' rwork size =',i4,'   iwork size =',i4/
     2  ' number of steps =',i5/
     3  ' number of f-s   =',i5/
     4  ' (excluding J-s) =',i5/
     5  ' number of J-s   =',i5/
     6  ' error overrun =',d10.2)
      if (miter .eq. 1 .or. miter .eq. 2)
     1   write (lout,185)nnz,ngp,nlu,nzl,nzu
 185  format(' number of nonzeros in J = ',i5/
     1  ' number of J index groups =',i5/
     2  ' number of LU decomp-s    =',i5/
     3  ' nonzeros in strict lower factor =',i5/
     4  ' nonzeros in strict upper factor =',i5)
      if (istate .lt. 0) go to 191
      if (miter .eq. 1 .or. miter .eq. 2)
     1   call ssout (neq, rwork(21), iwork, lout)
c Return for second run: minimal work array lengths, 1 output point.
      irun = irun + 1
      lrw = lenrw
      liw = leniw
      nout = 1
      go to 110
 191  continue
 192  continue
 193  continue
c
      write (lout,100)
      write (lout,200) nerr
 200  format(//'Number of errors encountered =',i3)
      stop
      end

      subroutine fdem (neq, t, y, ydot)
      integer neq,  i, igrid, j, l, m
      double precision t, y, ydot
      dimension y(neq), ydot(neq)
      data igrid/3/
      do 5 i = 1,neq
 5      ydot(i) = 0.0d0
      do 20 m = 1,igrid
        do 10 l = 1,igrid
          j = l + (m - 1)*igrid
          if (m .ne. 1) ydot(j-igrid) = ydot(j-igrid) + y(j)
          if (l .ne. 1) ydot(j-1) = ydot(j-1) + y(j)
          ydot(j) = ydot(j) - 4.0d0*y(j)
          if (l .ne. igrid) ydot(j+1) = ydot(j+1) + y(j)
 10       continue
 20     continue
      return
      end

      subroutine jdem (neq, t, y, j, ia, ja, pdj)
      integer neq, j, ia, ja,  igrid, l, m
      double precision t, y, pdj
      dimension y(neq), ia(*), ja(*), pdj(neq)
      data igrid/3/
      m = (j - 1)/igrid + 1
      l = j - (m - 1)*igrid
      pdj(j) = -4.0d0
      if (m .ne. 1) pdj(j-igrid) = 1.0d0
      if (l .ne. 1) pdj(j-1) = 1.0d0
      if (l .ne. igrid) pdj(j+1) = 1.0d0
      return
      end

      subroutine edit (y, iout, erm)
      integer iout,  i, neq
      double precision y, erm,   er, yex
      dimension y(*),yex(9,3)
      data neq /9/
      data yex /6.687279d-01, 9.901910d-01, 7.603061d-01,
     1   8.077979d-01, 1.170226e+00, 8.810605d-01, 5.013331d-01,
     2   7.201389d-01, 5.379644d-01, 1.340488d-01, 1.917157d-01,
     3   1.374034d-01, 1.007882d-01, 1.437868d-01, 1.028010d-01,
     4   3.844343d-02, 5.477593d-02, 3.911435d-02, 1.929166d-02,
     5   2.735444d-02, 1.939611d-02, 1.055981d-02, 1.496753d-02,
     6   1.060897d-02, 2.913689d-03, 4.128975d-03, 2.925977d-03/
      erm = 0.0d0
      do 10 i = 1,neq
        er = abs(y(i) - yex(i,iout))
 10     erm = max(erm,er)
      return
      end

      subroutine ssout (neq, iwk, iwork, lout)
      integer neq, iwk, iwork, lout
      integer i, i1, i2, ipian, ipjan, nnz
      dimension iwk(*), iwork(*)
      ipian = iwork(23)
      ipjan = iwork(24)
      nnz = iwork(19)
      i1 = ipian
      i2 = i1 + neq
      write (lout,10)(iwk(i),i=i1,i2)
 10   format(/' structure descriptor array ian ='/(20i4))
      i1 = ipjan
      i2 = i1 + nnz - 1
      write (lout,20)(iwk(i),i=i1,i2)
 20   format(/' structure descriptor array jan ='/(20i4))
      return
      end
