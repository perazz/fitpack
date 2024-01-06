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

    public :: FITPACK_SUCCESS_c
    public :: fitpack_message_c
    public :: curfit_c
    public :: percur_c
    public :: parcur_c
    public :: clocur_c
    public :: cocosp_c
    public :: concon_c

    contains

      ! Error test wrapper
      logical(FP_BOOL) function FITPACK_SUCCESS_c(ierr) result(success) bind(C,name="FITPACK_SUCCESS_c")
         integer(FP_FLAG), intent(in), value :: ierr
         success = FITPACK_SUCCESS(ierr)
      end function FITPACK_SUCCESS_c

      ! Flow control: on output flag present, return it; otherwise, halt on error
      subroutine fitpack_message_c(ierr,msg) bind(C,name='fitpack_message_c')
          integer(FP_FLAG), intent(in), value :: ierr
          character(len=1,kind=c_char), intent(inout) :: msg(*)
          character(:), allocatable :: str
          integer :: i

          str = FITPACK_MESSAGE(ierr)
          do i=1,len_trim(str)
             msg(i) = str(i:i)
          end do
          msg(len_trim(str)+1) = c_null_char
      end subroutine fitpack_message_c

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

      ! clocur interface
      subroutine clocur_c(iopt,ipar,idim,m,u,mx,x,w,k,s,nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier) bind(C,name='clocur_c')
          real   (FP_REAL), intent(in),    value :: s
          real   (FP_REAL), intent(inout)        :: fp
          integer(FP_SIZE), intent(in),    value :: iopt,ipar,idim,m,mx,k,nest,nc,lwrk
          integer(FP_SIZE), intent(inout)        :: n
          integer(FP_FLAG), intent(inout)        :: ier
          real   (FP_REAL), intent(in)           :: x(mx),w(m)
          real   (FP_REAL), intent(inout)        :: u(m),t(nest),c(nc),wrk(lwrk)
          integer(FP_SIZE), intent(inout)        :: iwrk(nest)
          call clocur(iopt,ipar,idim,m,u,mx,x,w,k,s,nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier)
      end subroutine clocur_c

      ! cocosp interface
      subroutine cocosp_c(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,wrk,lwrk,iwrk,kwrk,ier) bind(C,name='cocosp_c')
          real   (FP_REAL), intent(out)       :: sq
          integer(FP_SIZE), intent(in), value :: m,n,maxtr,maxbin,lwrk,kwrk
          integer(FP_FLAG), intent(out)       :: ier
          real   (FP_REAL), intent(in)        :: x(m),y(m),w(m),t(n)
          real   (FP_REAL), intent(inout)     :: e(n)
          real   (FP_REAL), intent(out)       :: c(n),sx(m)
          real   (FP_REAL), intent(inout)     :: wrk(lwrk)
          integer(FP_SIZE), intent(inout)     :: iwrk(kwrk)
          logical(FP_BOOL), intent(out)       :: bind(n)
          call cocosp(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,wrk,lwrk,iwrk,kwrk,ier)
      end subroutine cocosp_c

      ! concon interface
      subroutine concon_c(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx,bind,wrk,lwrk,iwrk,kwrk,ier) bind(C,name='concon_c')
          real   (FP_REAL), intent(in), value :: s
          real   (FP_REAL), intent(out)       :: sq
          integer(FP_SIZE), intent(in), value :: iopt,m,nest,maxtr,maxbin,lwrk,kwrk
          integer(FP_SIZE), intent(inout)     :: n
          integer(FP_FLAG), intent(out)       :: ier
          real   (FP_REAL), intent(in)        :: x(m),y(m),w(m)
          real   (FP_REAL), intent(inout)     :: v(m),t(nest),c(nest),sx(m),wrk(lwrk)
          integer(FP_SIZE), intent(inout)     :: iwrk(kwrk)
          logical(FP_BOOL), intent(inout)     :: bind(nest)
          call concon(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx,bind,wrk,lwrk,iwrk,kwrk,ier)
      end subroutine concon_c



end module fitpack_core_c
