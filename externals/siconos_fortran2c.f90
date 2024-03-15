! Siconos is a program dedicated to modeling, simulation and control
! of non smooth dynamical systems.
!
! Copyright 2024 INRIA.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
! http://www.apache.org/licenses/LICENSE-2.0
!
 ! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!

!> \file siconos_fortran2c.f90
!> A module used to wrap odepack functions with iso_c_binding stuff.

module siconos_fortran2c

  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: dlsodar2c, ql00012c, hem52c, n2qn12c
  
contains
    
  subroutine dlsodar2c(C_F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ISTATE, RWORK, LRW, IWORK, LIW,&
       C_JAC, JT, C_G, NG, JROOT) bind(c, name="lsodar")

    type(c_funptr), value :: C_F, C_G, C_JAC

    integer(c_int), intent(inout) :: NEQ
    real(c_double), intent(inout) :: Y(NEQ)
    real(c_double), intent(inout) :: T
    real(c_double), intent(inout) :: TOUT
    integer(c_int), intent(inout) :: ITOL
    integer(c_double), intent(inout), target :: RTOL, ATOL
    integer(c_size_t), intent(in) :: JT
    integer(c_int), intent(inout) :: ISTATE
    integer(c_int), intent(inout) :: LRW, LIW, NG
    real(c_double), intent(inout), target :: RWORK
    integer(c_int), intent(inout), target :: IWORK
    integer(c_int), intent(inout), target :: JROOT
    procedure(), pointer :: F_F, F_JAC, F_G
    integer :: ITASK = 1 ! for normal computation of output values of y at t = TOUT.
    integer :: IOPT = 0 ! means no optional inputs are being used
    call C_F_PROCPOINTER ( CPTR = C_F, FPTR = F_F )
    call C_F_PROCPOINTER ( CPTR = C_JAC, FPTR = F_JAC )
    call C_F_PROCPOINTER ( CPTR = C_G, FPTR = F_G )
    
    call DLSODAR(F_F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, F_JAC, JT, F_G, NG, JROOT)
    
  endsubroutine dlsodar2c
  
  subroutine n2qn12c (n, x, f, g, dxmin, df1, epsabs,mode, binf, bsup, iz, rz) bind(c, name="n2qn1")
    !! Note FP: since there are no comments in qnb.f, which is a full F77 file, it's quite hard to guess properly all the intent/types. This probably needs to
    !! be reviewed and check properly.
    integer(c_int), intent(in) :: n
    real(c_double), intent(inout) :: x(n)
    real(c_double), intent(inout) :: f
    real(c_double), intent(inout) :: g(n)
    real(c_double), intent(inout) :: dxmin(n)
    real(c_double), intent(inout) :: df1
    real(c_double), intent(inout) :: epsabs
    integer(c_int), intent(inout) :: mode
    real(c_double), intent(in) :: binf(n), bsup(n)
    integer(c_int), intent(inout), target :: iz(29)
    real(c_double), intent(inout) :: rz(45)

    integer :: imp = 3!!0
    integer :: io = 16
    integer :: iter = 500
    logical :: reverse = .true.
    integer :: nsim = 3*500
    call n2qn1(n, x, f, g, dxmin, df1, epsabs, imp, io, mode, iter, nsim, binf, bsup, iz, rz, reverse)

  end subroutine n2qn12c

  subroutine ql00012c(m,me,mmax,n,nmax,mnn,c,d,a,b,xl,xu,x,u,iout,ifail,iprint,war,lwar,iwar,liwar,eps) bind(c, name="ql0001")

    integer(c_long), intent(inout) :: m, me, mmax, n, nmax, mnn
    real(c_double), intent(inout), target :: c(nmax, n), d(n), a(mmax, n), b(mmax)
    real(c_double), intent(inout), target :: xl(n), xu(n), x(n), u(mnn)
    integer(c_long), intent(inout) :: iout, ifail, iprint, lwar, liwar
    real(c_double), intent(inout), target :: war(lwar)
    integer(c_long), intent(inout), target :: iwar(liwar)
    real(c_double), intent(inout) :: eps

    call ql0001(m,me,mmax,n,nmax,mnn,c,d,a,b,xl,xu,x,u,iout,ifail,iprint,war,lwar,iwar,liwar,eps)
    
  endsubroutine ql00012c
  
  subroutine hem52c(NQ,NV,NU,NL,C_FPROB,T,Q,V,U,A,RLAM,TEND,H,RTOL,ATOL,ITOL,C_SOLOUT,IOUT,WK,LWK,IWK,&
       LIWK,IDID) bind(c, name="hem5")

    integer(c_long), intent(inout) :: NQ, NV, NU, NL
    type(c_funptr), value :: C_FPROB, C_SOLOUT
    real(c_double), intent(inout), target :: Q(NQ), V(NV), U(NU), A(NV), RLAM(NL)
    real(c_double), intent(inout) :: T, TEND, H
    real(c_double), intent(inout), target :: RTOL(1), ATOL(1)
    integer(c_long), intent(inout) :: ITOL, IOUT
    integer(c_long), intent(inout) :: LWK, LIWK
    real(c_double), intent(inout), target :: WK(LWK)
    integer(c_long), intent(inout), target :: IWK(LIWK)
    integer(c_long), intent(inout) :: IDID
    
    procedure(), pointer :: F_FPROB, F_SOLOUT
    
    call C_F_PROCPOINTER ( CPTR = C_FPROB, FPTR = F_FPROB)
    call C_F_PROCPOINTER ( CPTR = C_SOLOUT, FPTR = F_SOLOUT)

    call hem5(NQ,NV,NU,NL,F_FPROB,T,Q,V,U,A,RLAM,TEND,H,RTOL,ATOL,ITOL,F_SOLOUT,IOUT,WK,LWK,IWK,LIWK,IDID)
    
  end subroutine hem52c


  
end module
