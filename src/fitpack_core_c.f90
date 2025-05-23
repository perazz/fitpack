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
    public :: splev_c
    public :: splder_c
    public :: spalde_c
    public :: curev_c
    public :: cualde_c
    public :: insert_c

    contains

      ! Error test wrapper
      logical(FP_BOOL) function FITPACK_SUCCESS_c(ierr) result(success) bind(C,name="FITPACK_SUCCESS_c")
         integer(FP_FLAG), intent(in), value :: ierr
         success = FITPACK_SUCCESS(ierr)
      end function FITPACK_SUCCESS_c

      ! Flow control: on output flag present, return it; otherwise, halt on error
      pure subroutine fitpack_message_c(ierr,msg) bind(C,name='fitpack_message_c')
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
      pure subroutine curfit_c(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier) bind(C,name="curfit_c")
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
      pure subroutine percur_c(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier) bind(C,name="percur_c")
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
      pure subroutine parcur_c(iopt,ipar,idim,m,u,mx,x,w,ub,ue,k,s,nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier) bind(C,name="parcur_c")
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
      pure subroutine clocur_c(iopt,ipar,idim,m,u,mx,x,w,k,s,nest,n,t,nc,c,fp,wrk,lwrk,iwrk,ier) bind(C,name='clocur_c')
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
      pure subroutine cocosp_c(m,x,y,w,n,t,e,maxtr,maxbin,c,sq,sx,bind,wrk,lwrk,iwrk,kwrk,ier) bind(C,name='cocosp_c')
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
      pure subroutine concon_c(iopt,m,x,y,w,v,s,nest,maxtr,maxbin,n,t,c,sq,sx,bind,wrk,lwrk,iwrk,kwrk,ier) bind(C,name='concon_c')
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

      ! splev interface
      pure subroutine splev_c(t,n,c,k,x,y,m,e,ier) bind(C,name='splev_c')
          integer(FP_SIZE), intent(in), value :: n, k, m
          real(FP_REAL),    intent(in)         :: t(n)
          real(FP_REAL),    intent(in)         :: c(n)
          real(FP_REAL),    intent(in)         :: x(m)
          real(FP_REAL),    intent(out)        :: y(m)
          integer(FP_FLAG), intent(in), value  :: e
          integer(FP_FLAG), intent(out)        :: ier
          call splev(t,n,c,k,x,y,m,e,ier)
      end subroutine splev_c
      
      ! splder interface
      pure subroutine splder_c(t,n,c,k,nu,x,y,m,e,wrk,ier) bind(C,name='splder_c')
          integer(FP_SIZE), intent(in), value  :: n,k,nu,m
          integer(FP_FLAG), intent(in), value  :: e
          integer(FP_FLAG), intent(out)        :: ier
          real(FP_REAL),    intent(in)         :: t(n),c(n),x(m)
          real(FP_REAL),    intent(out)        :: y(m)
          real(FP_REAL),    intent(inout)      :: wrk(n)      
          call splder(t,n,c,k,nu,x,y,m,e,wrk,ier)
      end subroutine splder_c
      
      ! spalde interface
      pure subroutine spalde_c(t,n,c,k1,x,d,ier) bind(C,name='spalde_c')
          integer(FP_SIZE), intent(in), value  :: n,k1
          integer(FP_FLAG), intent(out)        :: ier
          real(FP_REAL),    intent(in), value  :: x
          real(FP_REAL),    intent(in)         :: t(n),c(n)
          real(FP_REAL),    intent(out)        :: d(k1)
          call spalde(t,n,c,k1,x,d,ier)
      end subroutine spalde_c

      ! curev interface
      pure subroutine curev_c(idim,t,n,c,nc,k,u,m,x,mx,ier) bind(C,name='curev_c')
          integer(FP_SIZE), intent(in), value :: idim,n,nc,k,m,mx
          integer(FP_FLAG), intent(out)       :: ier
          real(FP_REAL),    intent(in)        :: t(n),c(nc),u(m)
          real(FP_REAL),    intent(out)       :: x(mx) ! mx -> assume 2d (idim,m)      
          call curev(idim,t,n,c,nc,k,u,m,x,mx,ier)
      end subroutine curev_c
      
      ! cualde interface
      pure subroutine cualde_c(idim,t,n,c,nc,k1,u,d,nd,ier) bind(C,name='cualde_c')
          integer(FP_SIZE), intent(in), value  :: idim,n,nc,k1,nd
          integer(FP_FLAG), intent(out)        :: ier
          real(FP_REAL),    intent(in), value  :: u
          real(FP_REAL),    intent(in)         :: t(n),c(nc)
          real(FP_REAL),    intent(out)        :: d(nd)
          call cualde(idim,t,n,c,nc,k1,u,d,nd,ier)
      end subroutine cualde_c
      
      ! insert interface
      pure subroutine insert_c(iopt,t,n,c,k,x,nest,ier) bind(C,name='insert_c')
          integer(FP_SIZE), intent(in), value  :: iopt,k,nest
          real(FP_REAL),    intent(inout)      :: t(nest),c(nest)          
          real(FP_REAL),    intent(in), value  :: x
          integer(FP_SIZE), intent(inout)      :: n
          integer(FP_FLAG), intent(out)        :: ier
          call insert_inplace(iopt,t,n,c,k,x,nest,ier)
      end subroutine insert_c
      
      ! splint interface
      real(FP_REAL) function splint_c(t,n,c,k,a,b,wrk) bind(C,name='splint_c')
          real(FP_REAL),    intent(in), value :: a,b
          integer(FP_SIZE), intent(in), value :: n,k
          real(FP_REAL),    intent(in)        :: t(n),c(n)
          real(FP_REAL),    intent(inout)     :: wrk(n)
          splint_c = splint(t,n,c,k,a,b,wrk)
      end function splint_c
      
      ! fourier coefficients
      pure subroutine fourco_c(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier) bind(C,name='fourco_c')
          integer(FP_SIZE), intent(in), value :: n,m
          integer(FP_FLAG), intent(out)       :: ier
          real(FP_REAL),    intent(in)        :: t(n),c(n),alfa(m)
          real(FP_REAL),    intent(inout)     :: wrk1(n),wrk2(n)
          real(FP_REAL),    intent(out)       :: ress(m),resc(m)
          call fourco(t,n,c,alfa,m,ress,resc,wrk1,wrk2,ier)
      end subroutine fourco_c

      ! spline roots
      pure subroutine sproot_c(t,n,c,zeros,mest,m,ier) bind(C,name='sproot_c')
          integer(FP_SIZE), intent(in), value  :: n,mest
          integer(FP_SIZE), intent(out)        :: m
          integer(FP_FLAG), intent(out)        :: ier
          real(FP_REAL),    intent(in)         :: t(n),c(n)
          real(FP_REAL),    intent(out)        :: zeros(mest)
          call sproot(t,n,c,zeros,mest,m,ier)
      end subroutine sproot_c

      ! surface fit
      pure subroutine surfit_c(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,nmax,eps,nx, &
                               tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier) bind(C,name='surfit_c')
          real(FP_REAL),    intent(in), value  :: xb,xe,yb,ye,s,eps
          real(FP_REAL),    intent(inout)      :: fp
          integer(FP_SIZE), intent(in), value  :: iopt,m,kx,ky,nxest,nyest,nmax,lwrk1,lwrk2,kwrk
          integer(FP_SIZE), intent(inout)      :: nx,ny
          integer(FP_FLAG), intent(out)        :: ier
          real(FP_REAL),    intent(in)         :: z(m),w(m)
          real(FP_REAL),    intent(inout)      :: x(m),y(m),tx(nmax),ty(nmax), &
                                                  c((nxest-kx-1)*(nyest-ky-1)),&
                                                  wrk1(lwrk1),wrk2(lwrk2)
          integer(FP_SIZE), intent(inout)      :: iwrk(kwrk)
          call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,nmax,eps,nx, &
                               tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
      end subroutine surfit_c

      ! regular grid fit
      pure subroutine regrid_c(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
                               nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier) bind(C,name='regrid_c')
          real(FP_REAL),    intent(in), value :: xb,xe,yb,ye,s
          real(FP_REAL),    intent(out)       :: fp
          integer(FP_SIZE), intent(in), value :: iopt,mx,my,kx,ky,nxest,nyest,lwrk,kwrk
          integer(FP_SIZE), intent(inout)     :: nx,ny
          integer(FP_FLAG), intent(out)       :: ier
          real(FP_REAL),    intent(in)        :: x(mx),y(my),z(mx*my)
          real(FP_REAL),    intent(inout)     :: tx(nxest),ty(nyest),c((nxest-kx-1)*(nyest-ky-1)),wrk(lwrk)
          integer(FP_SIZE), intent(inout)     :: iwrk(kwrk)
          call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
                               nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
      end subroutine regrid_c

      subroutine polar_c(iopt,m,x,y,z,w,rad,s,nuest,nvest,eps,nu,tu,nv,tv,u,v,c,&
                         fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier) bind(C,name='polar_c')
          real(FP_REAL),    intent(in), value :: s,eps
          real(FP_REAL),    intent(inout)     :: fp
          integer(FP_SIZE), intent(in), value :: m,nuest,nvest,lwrk1,lwrk2,kwrk
          integer(FP_SIZE), intent(out)       :: nu,nv
          integer(FP_FLAG), intent(out)       :: ier
          type(c_funptr),   intent(in), value :: rad
          real(FP_REAL),    intent(in)        :: x(m),y(m),z(m),w(m)
          real(FP_REAL),    intent(out)       :: u(m),v(m),tu(nuest),tv(nvest),c((nuest-4)*(nvest-4))
          real(FP_REAL),    intent(inout)     :: wrk1(lwrk1),wrk2(lwrk2)
          integer(FP_SIZE), intent(in)        :: iopt(3)
          integer(FP_SIZE), intent(inout)     :: iwrk(kwrk)
          procedure(fitpack_polar_boundary), pointer :: frad
          
          call c_f_procpointer(rad,frad)
          call polar(iopt,m,x,y,z,w,frad,s,nuest,nvest,eps,nu,tu,nv,tv,u,v,c,&
                              fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)

      end subroutine polar_c

      subroutine pogrid_c(iopt,ider,mu,u,mv,v,z,z0,r,s,nuest,nvest,nu,tu,nv,tv, &
                          c,fp,wrk,lwrk,iwrk,kwrk,ier) bind(C,name='pogrid_c')
          integer(FP_SIZE), intent(in)        :: iopt(3),ider(2)
          real(FP_REAL),    intent(in), value :: z0,r,s
          real(FP_REAL),    intent(inout)     :: fp
          integer(FP_SIZE), intent(in), value :: mu,mv,nuest,nvest,lwrk,kwrk
          integer(FP_FLAG), intent(out)       :: ier
          integer(FP_SIZE), intent(inout)     :: nu,nv
          integer(FP_SIZE), intent(inout)     :: iwrk(kwrk)
          real(FP_REAL), intent(in)           :: u(mu),v(mv),z(mu*mv)
          real(FP_REAL), intent(inout)        :: c((nuest-4)*(nvest-4)),tu(nuest),tv(nvest),wrk(lwrk)
          
          call pogrid(iopt,ider,mu,u,mv,v,z,z0,r,s,nuest,nvest,nu,tu,nv,tv, &
                      c,fp,wrk,lwrk,iwrk,kwrk,ier)
          
      end subroutine pogrid_c

      pure subroutine sphere_c(iopt,m,teta,phi,r,w,s,ntest,npest,eps,nt,tt,np,tp,c, &
                               fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier) bind(C,name='sphere_c')
          integer(FP_SIZE), intent(in), value :: iopt,m,ntest,npest,lwrk1,lwrk2,kwrk
          real(FP_REAL),    intent(in), value :: s,eps
          real(FP_REAL),    intent(inout)     :: fp
          integer(FP_SIZE), intent(inout)     :: nt,np
          integer(FP_FLAG), intent(out)       :: ier
          real(FP_REAL),    intent(in)        :: teta(m),phi(m),r(m),w(m)
          real(FP_REAL),    intent(inout)     :: tt(ntest),tp(npest),c((ntest-4)*(npest-4)), &
                                                 wrk1(lwrk1),wrk2(lwrk2)
          integer(FP_SIZE), intent(inout)     :: iwrk(kwrk)
          call sphere(iopt,m,teta,phi,r,w,s,ntest,npest,eps,nt,tt,np,tp,c, &
                      fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
      end subroutine sphere_c

      integer(FP_FLAG) function spgrid_c(iopt,ider,mu,u,mv,v,r,r0,r1,s,nuest,nvest,nu,tu,nv,tv, &
                             c,fp,wrk,lwrk,iwrk,kwrk) bind(C,name='spgrid_c') result(ier)

          real(FP_REAL),    intent(in), value :: r0,r1,s
          real(FP_REAL),    intent(inout)     :: fp
          integer(FP_SIZE), intent(in), value :: mu,mv,nuest,nvest,lwrk,kwrk
          integer(FP_SIZE), intent(inout)     :: nu,nv
          integer(FP_SIZE), intent(in)        :: iopt(3),ider(4)
          integer(FP_SIZE), intent(inout)     :: iwrk(kwrk)
          real(FP_REAL),    intent(in)        :: u(mu),v(mv),r(mu*mv)
          real(FP_REAL),    intent(inout)     :: c((nuest-4)*(nvest-4)),tu(nuest),tv(nvest),wrk(lwrk)                             
          
          call spgrid(iopt,ider,mu,u,mv,v,r,r0,r1,s,&
                      nuest,nvest,nu,tu,nv,tv,c,fp,wrk,lwrk,iwrk,kwrk,ier)
                             
      end function spgrid_c

      integer(FP_FLAG) function parsur_c(iopt,ipar,idim,mu,u,mv,v,f,s,nuest,nvest,nu,tu,nv,tv, &
                             c,fp,wrk,lwrk,iwrk,kwrk) result(ier) bind(C,name='parsur_c')

          real(FP_REAL),    intent(in), value  :: s
          real(FP_REAL),    intent(inout)      :: fp
          integer(FP_SIZE), intent(in), value  :: iopt,idim,mu,mv,nuest,nvest,lwrk,kwrk
          integer(FP_SIZE), intent(inout)      :: nu,nv
          real(FP_REAL),    intent(in)         :: u(mu),v(mv),f(mu*mv*idim)
          real(FP_REAL),    intent(inout)      :: tu(nuest),tv(nvest),c((nuest-4)*(nvest-4)*idim),wrk(lwrk)
          integer(FP_SIZE), intent(in)         :: ipar(2)
          integer(FP_SIZE), intent(inout)      :: iwrk(kwrk)
          
          call parsur(iopt,ipar,idim,mu,u,mv,v,f,s,nuest,nvest,nu,tu,nv,tv, &
                      c,fp,wrk,lwrk,iwrk,kwrk,ier)
          
      end function parsur_c

      integer(FP_FLAG) function bispeu_c(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk) &
                       result(ier) bind(C,name='bispeu_c')

          integer(FP_SIZE), intent(in), value  :: nx,ny,kx,ky,m,lwrk
          real(FP_REAL), intent(in)           :: tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(m),y(m)
          real(FP_REAL), intent(inout)        :: wrk(lwrk)
          real(FP_REAL), intent(out)          :: z(m)
          
          call bispeu(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk,ier)
          
      end function bispeu_c

      integer(FP_FLAG) function bispev_c(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,iwrk,kwrk) &
                       result(ier) bind(C,name='bispev_c')

          integer(FP_SIZE), intent(in), value  :: nx,ny,kx,ky,mx,my,lwrk,kwrk
          integer(FP_SIZE), intent(inout)     :: iwrk(kwrk)
          real(FP_REAL), intent(in)           :: tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my)
          real(FP_REAL), intent(out)          :: z(mx*my)
          real(FP_REAL), intent(inout)        :: wrk(lwrk)
          
          call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,iwrk,kwrk,ier)
          
      end function bispev_c

      integer(FP_FLAG) function parder_c(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk,iwrk,kwrk) &
                       result(ier) bind(C,name='parder_c')

          integer(FP_SIZE), intent(in), value  :: nx,ny,kx,ky,nux,nuy,mx,my,lwrk,kwrk
          integer(FP_SIZE), intent(inout)     :: iwrk(kwrk)
          real(FP_REAL), intent(in)           :: tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my)
          real(FP_REAL), intent(out)          :: z(mx*my)
          real(FP_REAL), intent(inout)        :: wrk(lwrk)
          
          call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk,iwrk,kwrk,ier)
          
      end function parder_c

end module fitpack_core_c
