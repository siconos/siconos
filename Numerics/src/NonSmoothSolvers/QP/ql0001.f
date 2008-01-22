C Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
      subroutine ql0001(m,me,mmax,n,nmax,mnn,c,d,a,b,xl,xu,
     1      x,u,iout,ifail,iprint,war,lwar,iwar,liwar,eps)
c
c**************************************************************************
c
c
c             solution of quadratic programming problems
c
c
c
c   ql0001 solves the quadratic programming problem
c
c   minimize        .5*x'*C*x + d'*x
c   subject to      A(j)*x  +  b(j)   =  0  ,  j=1,...,me
c                   A(j)*x  +  b(j)  >=  0  ,  j=me+1,...,m
c                   xl  <=  x  <=  xu
c   
c   Here C must be an n by n symmetric and positive matrix, d an n-dimensional
c   vector, A an m by n matrix and b an m-dimensional vector. The above
c   situation is indicated by iwar(1)=1. Alternatively, i.e. if iwar(1)=0,
c   the objective function matrix can also be provided in factorized form.
c   In this case, C is an upper triangular matrix.
c
c   The subroutine reorganizes some data so that the problem can be solved
c   by a modification of an algorithm proposed by Powell (1983).
c
c
c   usage:
c
c   call ql0001(m,me,mmax,n,nmax,mnn,c,d,a,b,xl,xu,x,u,iout,ifail,
c              iprint,war,lwar,iwar,liwar)
c
c
c   Definition of the parameters:
c
c   m :        total number of constraints.
c   me :       number of equality constraints.
c   mmax :     row dimension of a. mmax must be at least one and greater
c              than m.
c   n :        number of variables.
c   nmax :     row dimension of C. nmax must be greater or equal to n.
c   mnn :      must be equal to m + n + n.
c   c(nmax,nmax): objective function matrix which should be symmetric and
c              positive definite. If iwar(1) = 0, c is supposed to be the
c              choleskey-factor of another matrix, i.e. c is upper
c              triangular.
c   d(nmax) :  contains the constant vector of the objective function.
c   a(mmax,nmax): contains the data matrix of the linear constraints.
c   b(mmax) :  contains the constant data of the linear constraints.
c   xl(n),xu(n): contain the lower and upper bounds for the variables.
c   x(n) :     on return, x contains the optimal solution vector.
c   u(mnn) :   on return, u contains the lagrange multipliers. The first
c              m positions are reserved for the multipliers of the m
c              linear constraints and the subsequent ones for the 
c              multipliers of the lower and upper bounds. On successful
c              termination, all values of u with respect to inequalities 
c              and bounds should be greater or equal to zero.
c   iout :     integer indicating the desired output unit number, i.e.
c              all write-statements start with 'write(iout,... '.
c   ifail :    shows the termination reason.
c      ifail = 0 :   successful return.
c      ifail = 1 :   too many iterations (more than 40*(n+m)).
c      ifail = 2 :   accuracy insufficient to satisfy convergence
c                    criterion.
c      ifail = 5 :   length of a working array is too short.
c      ifail > 10 :  the constraints are inconsistent.
c   iprint :   output control.
c      iprint = 0 :  no output of ql0001.
c      iprint > 0 :  brief output in error cases.
c   war(lwar) : real working array. the length lwar should be grater than
c               3*nmax*nmax/2 + 10*nmax + 2*mmax.
c   iwar(liwar): integer working array. the length liwar should be at
c              least n.
c              if iwar(1)=1 initially, then the cholesky decomposition
c              which is required by the dual algorithm to get the first
c              unconstrained minimum of the objective function, is
c              performed internally. otherwise, i.e. if iwar(1)=0, then
c              it is assumed that the user provides the initial fac-
c              torization by himself and stores it in the upper trian-
c              gular part of the array c.
c
c   a named common-block  /cmache/eps   must be provided by the user,
c   where eps defines a guess for the underlying machine precision.
c
c
c   author (c): k. schittkowski,
c               mathematisches institut,
c               universitaet bayreuth,
c               95440 bayreuth,
c               germany, f.r.
c
c
c   version:    1.5  (june, 1991)
c
c
c*********************************************************************
c
c
      integer nmax,mmax,n,mnn,lwar,liwar
      dimension c(nmax,n),d(n),a(mmax,n),b(mmax),
     1      xl(n),xu(n),x(n),u(mnn),war(lwar),iwar(liwar)
      double precision c,d,a,b,x,xl,xu,u,war,diag,zero,
     1      eps,qpeps,ten
      integer m,me,iout,ifail,iprint,iwar,inw1,inw2,in,j,lw,mn,i,
     1      idiag,info,nact,maxit
      logical lql
c
c     intrinsic functions:  dsqrt
c
C      add eps in parameter Vincent Acary
c      common /cmache/eps
c
c     constant data
c
      lql=.false.
      if (iwar(1).eq.1) lql=.true.
      zero=0.0d+0
      ten=1.d+1
      maxit=40*(m+n)
      qpeps=eps
      inw1=1
      inw2=inw1+mmax
c
c     prepare problem data for execution
c
      if (m.le.0) goto 20
      in=inw1
      do 10 j=1,m
      war(in)=-b(j)
   10 in=in+1
   20 lw=3*nmax*nmax/2 + 10*nmax + m
      if ((inw2+lw).gt.lwar) goto 80
      if (liwar.lt.n) goto 81
      if (mnn.lt.m+n+n) goto 82
      mn=m+n
c
c     call of ql0002
c
      call ql0002(n,m,me,mmax,mn,mnn,nmax,lql,a,war(inw1),
     1    d,c,xl,xu,x,nact,iwar,maxit,qpeps,info,diag,
     2    war(inw2),lw)
c
c     test of matrix corrections
c
      ifail=0
      if (info.eq.1) goto 40
      if (info.eq.2) goto 90
      idiag=0
      if ((diag.gt.zero).and.(diag.lt.1000.0)) idiag=diag
      if ((iprint.gt.0).and.(idiag.gt.0))
     1   write(iout,1000) idiag
      if (info .lt. 0) goto  70
c
c     reorder multiplier
c
      do  50 j=1,mnn
   50 u(j)=zero
      in=inw2-1
      if (nact.eq.0) goto 30
      do  60 i=1,nact
      j=iwar(i)
      u(j)=war(in+i)
   60 continue
   30 continue
      return
c
c     error messages
c
   70 ifail=-info+10
      if ((iprint.gt.0).and.(nact.gt.0))
     1     write(iout,1100) -info,(iwar(i),i=1,nact)
      return
   80 ifail=5
      if (iprint .gt. 0) write(iout,1200)
      return
   81 ifail=5
      if (iprint .gt. 0) write(iout,1210)
      return
   82 ifail=5
      if (iprint .gt. 0) write(iout,1220)
      return
   40 ifail=1
      if (iprint.gt.0) write(iout,1300) maxit
      return
   90 ifail=2
      if (iprint.gt.0) write(iout,1400)
      return
c
c     format-instructions
c
 1000 format(/8x,28h***ql: matrix g was enlarged,i3,
     1        20h-times by unitmatrix)
 1100 format(/8x,18h***ql: constraint ,i5,
     1        19h not consistent to ,/,(10x,10i5))
 1200 format(/8x,21h***ql: lwar too small)
 1210 format(/8x,22h***ql: liwar too small)
 1220 format(/8x,20h***ql: mnn too small)
 1300 format(/8x,37h***ql: too many iterations (more than,i6,1h))
 1400 format(/8x,50h***ql: accuracy insufficient to attain convergence) 
      end
c
      subroutine ql0002(n,m,meq,mmax,mn,mnn,nmax,lql,a,b,grad,g,
     1      xl,xu,x,nact,iact,maxit,vsmall,info,diag,w,lw)
c
c**************************************************************************
c
c
c   this subroutine solves the quadratic programming problem 
c
c       minimize      grad'*x  +  0.5 * x*g*x
c       subject to    a(k)*x  =  b(k)   k=1,2,...,meq,
c                     a(k)*x >=  b(k)   k=meq+1,...,m,
c                     xl  <=  x  <=  xu
c
c   the quadratic programming method proceeds from an initial cholesky-
c   decomposition of the objective function matrix, to calculate the
c   uniquely determined minimizer of the unconstrained problem. 
c   successively all violated constraints are added to a working set 
c   and a minimizer of the objective function subject to all constraints
c   in this working set is computed. it is possible that constraints
c   have to leave the working set.
c
c
c   description of parameters:
c
c     n        : is the number of variables.
c     m        : total number of constraints.
c     meq      : number of equality contraints.
c     mmax     : row dimension of a, dimension of b. mmax must be at
c                least one and greater or equal to m.
c     mn       : must be equal m + n.
c     mnn      : must be equal m + n + n.
c     nmax     : row diemsion of g. must be at least n.
c     lql      : determines initial decomposition.
c        lql = .false.  : the upper triangular part of the matrix g
c                         contains initially the cholesky-factor of a suitable
c                         decomposition.
c        lql = .true.   : the initial cholesky-factorisation of g is to be
c                         performed by the algorithm.
c     a(mmax,nmax) : a is a matrix whose columns are the constraints normals.
c     b(mmax)  : contains the right hand sides of the constraints.
c     grad(n)  : contains the objective function vector grad.
c     g(nmax,n): contains the symmetric objective function matrix.
c     xl(n), xu(n): contain the lower and upper bounds for x.
c     x(n)     : vector of variables.
c     nact     : final number of active constraints.
c     iact(k) (k=1,2,...,nact): indices of the final active constraints.
c     info     : reason for the return from the subroutine.
c         info = 0 : calculation was terminated successfully.
c         info = 1 : maximum number of iterations attained.
c         info = 2 : accuracy is insufficient to maintain increasing
c                    function values.
c         info < 0 : the constraint with index abs(info) and the con-
c                    straints whose indices are iact(k), k=1,2,...,nact,
c                    are inconsistent.
c     maxit    : maximum number of iterations.
c     vsmall   : required accuracy to be achieved (e.g. in the order of the 
c                machine precision for small and well-conditioned problems).
c     diag     : on return diag is equal to the multiple of the unit matrix
c                that was added to g to achieve positive definiteness.
c     w(lw)    : the elements of w(.) are used for working space. the length
c                of w must not be less than (1.5*nmax*nmax + 10*nmax + m).
c                when info = 0 on return, the lagrange multipliers of the
c                final active constraints are held in w(k), k=1,2,...,nact.
c   the values of n, m, meq, mmax, mn, mnn and nmax and the elements of
c   a, b, grad and g are not altered.
c
c   the following integers are used to partition w:
c     the first n elements of w hold lagrange multiplier estimates.
c     w(iwz+i+(n-1)*j) holds the matrix element z(i,j).
c     w(iwr+i+0.5*j*(j-1)) holds the upper triangular matrix
c       element r(i,j). the subsequent n components of w may be
c       treated as an extra column of r(.,.).
c     w(iww-n+i) (i=1,2,...,n) are used for temporary storage.
c     w(iww+i) (i=1,2,...,n) are used for temporary storage.
c     w(iwd+i) (i=1,2,...,n) holds g(i,i) during the calculation.
c     w(iwx+i) (i=1,2,...,n) holds variables that will be used to
c       test that the iterations increase the objective function.
c     w(iwa+k) (k=1,2,...,m) usually holds the reciprocal of the
c       length of the k-th constraint, but its sign indicates
c       whether the constraint is active.
c
c   
c   author:    k. schittkowski,
c              mathematisches institut,
c              universitaet bayreuth,
c              8580 bayreuth,
c              germany, f.r.
c
c   author of original version:
c              m.j.d. powell, damtp,
c              university of cambridge, silver street
c              cambridge,
c              england
c
c
c   reference: m.j.d. powell: zqpcvx, a fortran subroutine for convex
c              programming, report damtp/1983/na17, university of
c              cambridge, england, 1983.
c
c
c   version :  2.0 (march, 1987)
c
c
c*************************************************************************
c
      integer mmax,nmax,n,lw,nflag,iwwn
      dimension a(mmax,n),b(mmax),grad(n),g(nmax,n),x(n),iact(n),
     1      w(lw),xl(n),xu(n)
      integer m,meq,mn,mnn,nact,iact,info,maxit
      double precision cvmax,diag,diagr,fdiff,fdiffa,ga,gb,parinc,parnew
     1      ,ratio,res,step,sum,sumx,sumy,suma,sumb,sumc,temp,tempa,
     2       vsmall,xmag,xmagr,zero,one,two,onha,vfact
      double precision a,b,g,grad,w,x,xl,xu
c
c   intrinsic functions:   dmax1,dsqrt,dabs,dmin1
c
      integer iwz,iwr,iww,iwd,iwa,ifinc,kfinc,k,i,ia,id,ii,ir,ira,
     1     irb,j,nm,iz,iza,iterc,itref,jfinc,iflag,iws,is,k1,iw,kk,ip,
     2     ipp,il,iu,ju,kflag,lflag,jflag,kdrop,nu,mflag,knext,ix,iwx,
     3     iwy,iy,jl
      logical lql,lower
c
c   initial addresses
c
      iwz=nmax
      iwr=iwz+nmax*nmax
      iww=iwr+(nmax*(nmax+3))/2
      iwd=iww+nmax
      iwx=iwd+nmax
      iwa=iwx+nmax
c
c     set some constants.
c
      zero=0.d+0
      one=1.d+0
      two=2.d+0
      onha=1.5d+0
      vfact=1.d+0
c
c     set some parameters.
c     number less than vsmall are assumed to be negligible.
c     the multiple of i that is added to g is at most diagr times
c       the least multiple of i that gives positive definiteness.
c     x is re-initialised if its magnitude is reduced by the
c       factor xmagr.
c     a check is made for an increase in f every ifinc iterations,
c       after kfinc iterations are completed.
c
      diagr=two
      diag=zero
      xmagr=1.0d-2
      ifinc=3
      kfinc=max0(10,n)
c
c     find the reciprocals of the lengths of the constraint normals.
c     return if a constraint is infeasible due to a zero normal.
c
      nact=0
      if (m .le. 0) goto 45
      do 40 k=1,m
      sum=zero
      do 10 i=1,n
   10 sum=sum+a(k,i)**2
      if (sum .gt. zero) goto 20
      if (b(k) .eq. zero) goto 30
      info=-k
      if (k .le. meq) goto 730
      if (b(k)) 30,30,730
   20 sum=one/dsqrt(sum)
   30 ia=iwa+k
   40 w(ia)=sum
   45 do 50 k=1,n
      ia=iwa+m+k
   50 w(ia)=one
c
c     if necessary increase the diagonal elements of g.
c
      if (.not. lql) goto 165
      do 60 i=1,n
      id=iwd+i
      w(id)=g(i,i)
      diag=dmax1(diag,vsmall-w(id))
      if (i .eq. n) goto 60
      ii=i+1
      do 55 j=ii,n
      ga=-dmin1(w(id),g(j,j))
      gb=dabs(w(id)-g(j,j))+dabs(g(i,j))
      if (gb .gt. zero) ga=ga+g(i,j)**2/gb
   55 diag=dmax1(diag,ga)
   60 continue
      if (diag .le. zero) goto 90
   70 diag=diagr*diag
      do 80 i=1,n
      id=iwd+i
   80 g(i,i)=diag+w(id)
c
c     form the cholesky factorisation of g. the transpose
c     of the factor will be placed in the r-partition of w.
c
   90 ir=iwr
      do 130 j=1,n
      ira=iwr
      irb=ir+1
      do 120 i=1,j
      temp=g(i,j)
      if (i .eq. 1) goto 110
      do 100 k=irb,ir
      ira=ira+1
  100 temp=temp-w(k)*w(ira)
  110 ir=ir+1
      ira=ira+1
      if (i .lt. j) w(ir)=temp/w(ira)
  120 continue
      if (temp .lt. vsmall) goto 140
  130 w(ir)=dsqrt(temp)
      goto 170
c
c     increase further the diagonal element of g.
c
  140 w(j)=one
      sumx=one
      k=j
  150 sum=zero
      ira=ir-1
      do 160 i=k,j
      sum=sum-w(ira)*w(i)
  160 ira=ira+i
      ir=ir-k
      k=k-1
      w(k)=sum/w(ir)
      sumx=sumx+w(k)**2
      if (k .ge. 2) goto 150
      diag=diag+vsmall-temp/sumx
      goto 70
c
c     store the cholesky factorisation in the r-partition
c     of w.
c
  165 ir=iwr
      do 166 i=1,n
      do 166 j=1,i
      ir=ir+1
  166 w(ir)=g(j,i)
c
c     set z the inverse of the matrix in r.
c
  170 nm=n-1
      do 220 i=1,n
      iz=iwz+i
      if (i .eq. 1) goto 190
      do 180 j=2,i
      w(iz)=zero
  180 iz=iz+n
  190 ir=iwr+(i+i*i)/2
      w(iz)=one/w(ir)
      if (i .eq. n) goto 220
      iza=iz
      do 210 j=i,nm
      ir=ir+i
      sum=zero
      do 200 k=iza,iz,n
      sum=sum+w(k)*w(ir)
  200 ir=ir+1
      iz=iz+n
  210 w(iz)=-sum/w(ir)
  220 continue
c
c     set the initial values of some variables.
c     iterc counts the number of iterations.
c     itref is set to one when iterative refinement is required.
c     jfinc indicates when to test for an increase in f.
c
      iterc=1
      itref=0
      jfinc=-kfinc
c
c     set x to zero and set the corresponding residuals of the
c     kuhn-tucker conditions.
c
  230 iflag=1
      iws=iww-n
      do 240 i=1,n
      x(i)=zero
      iw=iww+i
      w(iw)=grad(i)
      if (i .gt. nact) goto 240
      w(i)=zero
      is=iws+i
      k=iact(i)
      if (k .le. m) goto 235
      if (k .gt. mn) goto 234
      k1=k-m
      w(is)=xl(k1)
      goto 240
  234 k1=k-mn
      w(is)=-xu(k1)
      goto 240
  235 w(is)=b(k)
  240 continue
      xmag=zero
      vfact=1.d+0
      if (nact) 340,340,280
c
c     set the residuals of the kuhn-tucker conditions for general x.
c
  250 iflag=2
      iws=iww-n
      do 260 i=1,n
      iw=iww+i
      w(iw)=grad(i)
      if (lql) goto 259
      id=iwd+i
      w(id)=zero
      do 251 j=i,n
  251 w(id)=w(id)+g(i,j)*x(j)
      do 252 j=1,i
      id=iwd+j
  252 w(iw)=w(iw)+g(j,i)*w(id)
      goto 260
  259 do 261 j=1,n
  261 w(iw)=w(iw)+g(i,j)*x(j)
  260 continue
      if (nact .eq. 0) goto 340
      do 270 k=1,nact
      kk=iact(k)
      is=iws+k
      if (kk .gt. m) goto 265
      w(is)=b(kk)
      do 264 i=1,n
      iw=iww+i
      w(iw)=w(iw)-w(k)*a(kk,i)
  264 w(is)=w(is)-x(i)*a(kk,i)
      goto 270
  265 if (kk .gt. mn) goto 266
      k1=kk-m
      iw=iww+k1
      w(iw)=w(iw)-w(k)
      w(is)=xl(k1)-x(k1)
      goto 270
  266 k1=kk-mn
      iw=iww+k1
      w(iw)=w(iw)+w(k)
      w(is)=-xu(k1)+x(k1)
  270 continue
c
c     pre-multiply the vector in the s-partition of w by the
c     invers of r transpose.
c
  280 ir=iwr
      ip=iww+1
      ipp=iww+n
      il=iws+1
      iu=iws+nact
      do 310 i=il,iu
      sum=zero
      if (i .eq. il) goto 300
      ju=i-1
      do 290 j=il,ju
      ir=ir+1
  290 sum=sum+w(ir)*w(j)
  300 ir=ir+1
  310 w(i)=(w(i)-sum)/w(ir)
c
c     shift x to satisfy the active constraints and make the
c     corresponding change to the gradient residuals.
c
      do 330 i=1,n
      iz=iwz+i
      sum=zero
      do 320 j=il,iu
      sum=sum+w(j)*w(iz)
  320 iz=iz+n
      x(i)=x(i)+sum
      if (lql) goto 329
      id=iwd+i
      w(id)=zero
      do 321 j=i,n
  321 w(id)=w(id)+g(i,j)*sum
      iw=iww+i
      do 322 j=1,i
      id=iwd+j
  322 w(iw)=w(iw)+g(j,i)*w(id)
      goto 330
  329 do 331 j=1,n
      iw=iww+j
  331 w(iw)=w(iw)+sum*g(i,j)
  330 continue
c
c     form the scalar product of the current gradient residuals
c     with each column of z.
c
  340 kflag=1
      goto 930
  350 if (nact .eq. n) goto 380
c
c     shift x so that it satisfies the remaining kuhn-tucker
c     conditions.
c
      il=iws+nact+1
      iza=iwz+nact*n
      do 370 i=1,n
      sum=zero
      iz=iza+i
      do 360 j=il,iww
      sum=sum+w(iz)*w(j)
  360 iz=iz+n
  370 x(i)=x(i)-sum
      info=0
      if (nact .eq. 0) goto 410
c
c     update the lagrange multipliers.
c
  380 lflag=3
      goto 740
  390 do 400 k=1,nact
      iw=iww+k
  400 w(k)=w(k)+w(iw)
c
c     revise the values of xmag.
c     branch if iterative refinement is required.
c
  410 jflag=1
      goto 910
  420 if (iflag .eq. itref) goto 250
c
c     delete a constraint if a lagrange multiplier of an
c     inequality constraint is negative.
c
      kdrop=0
      goto 440
  430 kdrop=kdrop+1
      if (w(kdrop) .ge. zero) goto 440
      if (iact(kdrop) .le. meq) goto 440
      nu=nact
      mflag=1
      goto 800
  440 if (kdrop .lt. nact) goto 430
c
c     seek the greateast normalised constraint violation, disregarding
c     any that may be due to computer rounding errors.
c
  450 cvmax=zero
      if (m .le. 0) goto 481
      do 480 k=1,m
      ia=iwa+k
      if (w(ia) .le. zero) goto 480
      sum=-b(k)
      do 460 i=1,n
  460 sum=sum+x(i)*a(k,i)
      sumx=-sum*w(ia)
      if (k .le. meq) sumx=dabs(sumx)
      if (sumx .le. cvmax) goto 480
      temp=dabs(b(k))
      do 470 i=1,n
  470 temp=temp+dabs(x(i)*a(k,i))
      tempa=temp+dabs(sum)
      if (tempa .le. temp) goto 480
      temp=temp+onha*dabs(sum)
      if (temp .le. tempa) goto 480
      cvmax=sumx
      res=sum
      knext=k
  480 continue
  481 do 485 k=1,n
      lower=.true.
      ia=iwa+m+k
      if (w(ia) .le. zero) goto 485
      sum=xl(k)-x(k)
      if (sum) 482,485,483
  482 sum=x(k)-xu(k)
      lower=.false.
  483 if (sum .le. cvmax) goto 485
      cvmax=sum
      res=-sum
      knext=k+m
      if (lower) goto 485
      knext=k+mn
  485 continue
c
c     test for convergence
c
      info=0
      if (cvmax .le. vsmall) goto 700
c
c     return if, due to rounding errors, the actual change in
c     x may not increase the objective function
c
      jfinc=jfinc+1
      if (jfinc .eq. 0) goto 510
      if (jfinc .ne. ifinc) goto 530
      fdiff=zero
      fdiffa=zero
      do 500 i=1,n
      sum=two*grad(i)
      sumx=dabs(sum)
      if (lql) goto 489
      id=iwd+i
      w(id)=zero
      do 486 j=i,n
      ix=iwx+j
  486 w(id)=w(id)+g(i,j)*(w(ix)+x(j))
      do 487 j=1,i
      id=iwd+j
      temp=g(j,i)*w(id)
      sum=sum+temp
  487 sumx=sumx+dabs(temp)
      goto 495
  489 do 490 j=1,n
      ix=iwx+j
      temp=g(i,j)*(w(ix)+x(j))
      sum=sum+temp
  490 sumx=sumx+dabs(temp)
  495 ix=iwx+i
      fdiff=fdiff+sum*(x(i)-w(ix))
  500 fdiffa=fdiffa+sumx*dabs(x(i)-w(ix))
      info=2
      sum=fdiffa+fdiff
      if (sum .le. fdiffa) goto 700
      temp=fdiffa+onha*fdiff
      if (temp .le. sum) goto 700
      jfinc=0
      info=0
  510 do 520 i=1,n
      ix=iwx+i
  520 w(ix)=x(i)
c
c     form the scalar product of the new constraint normal with each
c     column of z. parnew will become the lagrange multiplier of
c     the new constraint.
c
  530 iterc=iterc+1
      if (iterc.le.maxit) goto 531
      info=1
      goto 710
  531 continue
      iws=iwr+(nact+nact*nact)/2
      if (knext .gt. m) goto 541
      do 540 i=1,n
      iw=iww+i
  540 w(iw)=a(knext,i)
      goto 549
  541 do 542 i=1,n
      iw=iww+i
  542 w(iw)=zero
      k1=knext-m
      if (k1 .gt. n) goto 545
      iw=iww+k1
      w(iw)=one
      iz=iwz+k1
      do 543 i=1,n
      is=iws+i
      w(is)=w(iz)
  543 iz=iz+n
      goto 550
  545 k1=knext-mn
      iw=iww+k1
      w(iw)=-one
      iz=iwz+k1
      do 546 i=1,n
      is=iws+i
      w(is)=-w(iz)
  546 iz=iz+n
      goto 550
  549 kflag=2
      goto 930
  550 parnew=zero
c
c     apply givens rotations to make the last (n-nact-2) scalar
c     products equal to zero.
c
      if (nact .eq. n) goto 570
      nu=n
      nflag=1
      goto 860
c
c     branch if there is no need to delete a constraint.
c
  560 is=iws+nact
      if (nact .eq. 0) goto 640
      suma=zero
      sumb=zero
      sumc=zero
      iz=iwz+nact*n
      do 563 i=1,n
      iz=iz+1
      iw=iww+i
      suma=suma+w(iw)*w(iz)
      sumb=sumb+dabs(w(iw)*w(iz))
  563 sumc=sumc+w(iz)**2
      temp=sumb+.1d+0*dabs(suma)
      tempa=sumb+.2d+0*dabs(suma)
      if (temp .le. sumb) goto 570
      if (tempa .le. temp) goto 570
      if (sumb .gt. vsmall) goto 5
      goto 570
    5 sumc=dsqrt(sumc)
      ia=iwa+knext
      if (knext .le. m) sumc=sumc/w(ia)
      temp=sumc+.1d+0*dabs(suma)
      tempa=sumc+.2d+0*dabs(suma)
      if (temp .le. sumc) goto 567
      if (tempa .le. temp) goto 567
      goto 640
c
c     calculate the multipliers for the new constraint normal
c     expressed in terms of the active constraint normals.
c     then work out which contraint to drop.
c
  567 lflag=4
      goto 740
  570 lflag=1
      goto 740
c
c     complete the test for linearly dependent constraints.
c
  571 if (knext .gt. m) goto 574
      do 573 i=1,n
      suma=a(knext,i)
      sumb=dabs(suma)
      if (nact.eq.0) goto 581
      do 572 k=1,nact
      kk=iact(k)
      if (kk.le.m) goto 568
      kk=kk-m
      temp=zero
      if (kk.eq.i) temp=w(iww+kk)
      kk=kk-n
      if (kk.eq.i) temp=-w(iww+kk)
      goto 569
  568 continue
      iw=iww+k
      temp=w(iw)*a(kk,i)
  569 continue
      suma=suma-temp
  572 sumb=sumb+dabs(temp)
  581 if (suma .le. vsmall) goto 573
      temp=sumb+.1d+0*dabs(suma)
      tempa=sumb+.2d+0*dabs(suma)
      if (temp .le. sumb) goto 573
      if (tempa .le. temp) goto 573
      goto 630
  573 continue
      lflag=1
      goto 775
  574 k1=knext-m
      if (k1 .gt. n) k1=k1-n
      do 578 i=1,n
      suma=zero
      if (i .ne. k1) goto 575
      suma=one
      if (knext .gt. mn) suma=-one
  575 sumb=dabs(suma)
      if (nact.eq.0) goto 582
      do 577 k=1,nact
      kk=iact(k)
      if (kk .le. m) goto 579
      kk=kk-m
      temp=zero
      if (kk.eq.i) temp=w(iww+kk)
      kk=kk-n
      if (kk.eq.i) temp=-w(iww+kk)
      goto 576
  579 iw=iww+k
      temp=w(iw)*a(kk,i)
  576 suma=suma-temp
  577 sumb=sumb+dabs(temp)
  582 temp=sumb+.1d+0*dabs(suma)
      tempa=sumb+.2d+0*dabs(suma)
      if (temp .le. sumb) goto 578
      if (tempa .le. temp) goto 578
      goto 630
  578 continue
      lflag=1
      goto 775
c
c     branch if the contraints are inconsistent.
c
  580 info=-knext
      if (kdrop .eq. 0) goto 700
      parinc=ratio
      parnew=parinc
c
c     revise the lagrange multipliers of the active constraints.
c
  590 if (nact.eq.0) goto 601
      do 600 k=1,nact
      iw=iww+k
      w(k)=w(k)-parinc*w(iw)
      if (iact(k) .gt. meq) w(k)=dmax1(zero,w(k))
  600 continue
  601 if (kdrop .eq. 0) goto 680
c
c     delete the constraint to be dropped.
c     shift the vector of scalar products.
c     then, if appropriate, make one more scalar product zero.
c
      nu=nact+1
      mflag=2
      goto 800
  610 iws=iws-nact-1
      nu=min0(n,nu)
      do 620 i=1,nu
      is=iws+i
      j=is+nact
  620 w(is)=w(j+1)
      nflag=2
      goto 860
c
c     calculate the step to the violated constraint.
c
  630 is=iws+nact
  640 sumy=w(is+1)
      step=-res/sumy
      parinc=step/sumy
      if (nact .eq. 0) goto 660
c
c     calculate the changes to the lagrange multipliers, and reduce
c     the step along the new search direction if necessary.
c
      lflag=2
      goto 740
  650 if (kdrop .eq. 0) goto 660
      temp=one-ratio/parinc
      if (temp .le. zero) kdrop=0
      if (kdrop .eq. 0) goto 660
      step=ratio*sumy
      parinc=ratio
      res=temp*res
c
c     update x and the lagrange multipiers.
c     drop a constraint if the full step is not taken.
c
  660 iwy=iwz+nact*n
      do 670 i=1,n
      iy=iwy+i
  670 x(i)=x(i)+step*w(iy)
      parnew=parnew+parinc
      if (nact .ge. 1) goto 590
c
c     add the new constraint to the active set.
c
  680 nact=nact+1
      w(nact)=parnew
      iact(nact)=knext
      ia=iwa+knext
      if (knext .gt. mn) ia=ia-n
      w(ia)=-w(ia)
c
c     estimate the magnitude of x. then begin a new iteration,
c     re-initilising x if this magnitude is small.
c
      jflag=2
      goto 910
  690 if (sum .lt. (xmagr*xmag)) goto 230
      if (itref) 450,450,250
c
c     initiate iterative refinement if it has not yet been used,
c     or return after restoring the diagonal elements of g.
c
  700 if (iterc .eq. 0) goto 710
      itref=itref+1
      jfinc=-1
      if (itref .eq. 1) goto 250
  710 if (.not. lql) return
      do 720 i=1,n
      id=iwd+i
  720 g(i,i)=w(id)
  730 return
c
c
c     the remaining instructions are used as subroutines.
c
c
c********************************************************************
c
c
c     calculate the lagrange multipliers by pre-multiplying the
c     vector in the s-partition of w by the inverse of r.
c
  740 ir=iwr+(nact+nact*nact)/2
      i=nact
      sum=zero
      goto 770
  750 ira=ir-1
      sum=zero
      if (nact.eq.0) goto 761
      do 760 j=i,nact
      iw=iww+j
      sum=sum+w(ira)*w(iw)
  760 ira=ira+j
  761 ir=ir-i
      i=i-1
  770 iw=iww+i
      is=iws+i
      w(iw)=(w(is)-sum)/w(ir)
      if (i .gt. 1) goto 750
      if (lflag .eq. 3) goto 390
      if (lflag .eq. 4) goto 571
c
c     calculate the next constraint to drop.
c
  775 ip=iww+1
      ipp=iww+nact
      kdrop=0
      if (nact.eq.0) goto 791
      do 790 k=1,nact
      if (iact(k) .le. meq) goto 790
      iw=iww+k
      if ((res*w(iw)) .ge. zero) goto 790
      temp=w(k)/w(iw)
      if (kdrop .eq. 0) goto 780
      if (dabs(temp) .ge. dabs(ratio)) goto 790
  780 kdrop=k
      ratio=temp
  790 continue
  791 goto (580,650), lflag
c
c
c********************************************************************
c
c
c     drop the constraint in position kdrop in the active set.
c
  800 ia=iwa+iact(kdrop)
      if (iact(kdrop) .gt. mn) ia=ia-n
      w(ia)=-w(ia)
      if (kdrop .eq. nact) goto 850
c
c     set some indices and calculate the elements of the next
c     givens rotation.
c
      iz=iwz+kdrop*n
      ir=iwr+(kdrop+kdrop*kdrop)/2
  810 ira=ir
      ir=ir+kdrop+1
      temp=dmax1(dabs(w(ir-1)),dabs(w(ir)))
      sum=temp*dsqrt((w(ir-1)/temp)**2+(w(ir)/temp)**2)
      ga=w(ir-1)/sum
      gb=w(ir)/sum
c
c     exchange the columns of r.
c
      do 820 i=1,kdrop
      ira=ira+1
      j=ira-kdrop
      temp=w(ira)
      w(ira)=w(j)
  820 w(j)=temp
      w(ir)=zero
c
c     apply the rotation to the rows of r.
c
      w(j)=sum
      kdrop=kdrop+1
      do 830 i=kdrop,nu
      temp=ga*w(ira)+gb*w(ira+1)
      w(ira+1)=ga*w(ira+1)-gb*w(ira)
      w(ira)=temp
  830 ira=ira+i
c
c     apply the rotation to the columns of z.
c
      do 840 i=1,n
      iz=iz+1
      j=iz-n
      temp=ga*w(j)+gb*w(iz)
      w(iz)=ga*w(iz)-gb*w(j)
  840 w(j)=temp
c
c     revise iact and the lagrange multipliers.
c
      iact(kdrop-1)=iact(kdrop)
      w(kdrop-1)=w(kdrop)
      if (kdrop .lt. nact) goto 810
  850 nact=nact-1
      goto (250,610), mflag
c
c
c********************************************************************
c
c
c     apply givens rotation to reduce some of the scalar
c     products in the s-partition of w to zero.
c
  860 iz=iwz+nu*n
  870 iz=iz-n
  880 is=iws+nu
      nu=nu-1
      if (nu .eq. nact) goto 900
      if (w(is) .eq. zero) goto 870
      temp=dmax1(dabs(w(is-1)),dabs(w(is)))
      sum=temp*dsqrt((w(is-1)/temp)**2+(w(is)/temp)**2)
      ga=w(is-1)/sum
      gb=w(is)/sum
      w(is-1)=sum
      do 890 i=1,n
      k=iz+n
      temp=ga*w(iz)+gb*w(k)
      w(k)=ga*w(k)-gb*w(iz)
      w(iz)=temp
  890 iz=iz-1
      goto 880
  900 goto (560,630), nflag
c
c
c********************************************************************
c
c
c     calculate the magnitude of x an revise xmag.
c
  910 sum=zero
      do 920 i=1,n
      sum=sum+dabs(x(i))*vfact*(dabs(grad(i))+dabs(g(i,i)*x(i)))
      if (lql) goto 920
      if (sum .lt. 1.d-30) goto 920
      vfact=1.d-10*vfact
      sum=1.d-10*sum
      xmag=1.d-10*xmag
  920 continue
  925 xmag=dmax1(xmag,sum)
      goto (420,690), jflag
c
c
c********************************************************************
c
c
c     pre-multiply the vector in the w-partition of w by z transpose.
c
  930 jl=iww+1
      iz=iwz
      do 940 i=1,n
      is=iws+i
      w(is)=zero
      iwwn=iww+n
      do 940 j=jl,iwwn
      iz=iz+1
  940 w(is)=w(is)+w(iz)*w(j)
      goto (350,550), kflag
      return
      end
