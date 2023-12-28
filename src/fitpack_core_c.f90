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

      ! Error test wrapper
      logical(FP_BOOL) function FITPACK_SUCCESS_c(ierr) result(success) bind(C,name="FITPACK_SUCCESS_c")
         integer(FP_FLAG), intent(in), value :: ierr
         success = FITPACK_SUCCESS(ierr)
      end function FITPACK_SUCCESS_c

      ! curfit interface
      subroutine curfit_c(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier) bind(C,name="curfit_c")
         real   (FP_REAL), intent(in),    value  :: xb,xe,s
         real   (FP_REAL), intent(inout)         :: fp
         integer(FP_SIZE), intent(in),    value  :: iopt,m,k,nest,lwrk
         integer(FP_FLAG), intent(out)           :: ier
         integer(FP_SIZE), intent(inout)         :: n
         real   (FP_REAL), intent(in)            :: x(m),y(m),w(m)
         real   (FP_REAL), intent(inout)         :: t(nest),c(nest),wrk(lwrk)
         integer(FP_SIZE), intent(inout)         :: iwrk(nest)
         call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
      end subroutine curfit_c

      ! percur interface
      subroutine percur_c(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier) bind(C,name="percur_c")
         real   (FP_REAL), intent(in),    value :: s
         real   (FP_REAL), intent(inout)        :: fp
         integer(FP_SIZE), intent(inout)        :: n
         integer(FP_FLAG), intent(inout)        :: ier
         integer(FP_SIZE), intent(in),    value :: iopt,m,k,nest,lwrk
         real   (FP_REAL), intent(in)           :: x(m),y(m),w(m)
         real   (FP_REAL), intent(inout)        :: t(nest),c(nest),wrk(lwrk)
         integer(FP_SIZE), intent(inout)        :: iwrk(nest)
         call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
      end subroutine percur_c

      ! parcur interface
      subroutine parcur_c(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier) bind(C,name="parcur_c")
          real   (FP_REAL), intent(inout)       :: ub,ue,s
          real   (FP_REAL), intent(out)         :: fp
          integer(FP_SIZE), intent(in),   value :: iopt,ipar,idim,m,mx,k,nest,lwrk,nc
          integer(FP_SIZE), intent(inout)       :: n
          integer(FP_FLAG), intent(out)         :: ier
          real   (FP_REAL), intent(in)          :: x(idim,m)
          real   (FP_REAL), intent(inout)       :: u(m),w(m),t(nest),c(nc),wrk(lwrk)
          integer(FP_SIZE), intent(inout)       :: iwrk(nest)
          call parcur(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier)
      end subroutine parcur_c


end module fitpack_core_c
