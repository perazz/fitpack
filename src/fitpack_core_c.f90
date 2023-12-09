! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   (C) Federico Perini, 12/09/2023
!   Based on the netlib library by Paul Dierckx
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************
module fitpack_core_c
    use fitpack_core
    use iso_c_binding
    implicit none
    private

    contains

      ! curfit interface
      subroutine curfit_c(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier) bind(C,name="curfit_c")
         real(FP_REAL),    intent(in),    value  :: xb,xe,s
         real(FP_REAL),    intent(inout)         :: fp
         integer(FP_SIZE), intent(in),    value  :: iopt,m,k,nest,lwrk
         integer(FP_FLAG), intent(out)           :: ier
         integer(FP_SIZE), intent(inout)         :: n
         real(FP_REAL),    intent(in)            :: x(m),y(m),w(m)
         real(FP_REAL),    intent(inout)         :: t(nest),c(nest),wrk(lwrk)
         integer(FP_SIZE), intent(inout)         :: iwrk(nest)
         call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
      end subroutine curfit_c


end module fitpack_core_c
