! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   Refactored by Federico Perini, 10/6/2022
!   Based on the netlib library by Paul Dierckx
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************
module fitpack_fitters
    use fitpack_core
    implicit none
    private

    public :: fitpack_fitter

    !> Abstract base type for all fitpack OOP fitters.
    !> Contains the shared fitting state fields and communication helpers.
    type, abstract :: fitpack_fitter

        !> Fitting state flag
        integer(FP_FLAG) :: iopt = IOPT_NEW_SMOOTHING

        !> Smoothing parameter
        real(FP_REAL) :: smoothing = 1000.0_FP_REAL

        !> Weighted sum of squared residuals
        real(FP_REAL) :: fp = zero

        !> B-spline coefficients
        real(FP_REAL), allocatable :: c(:)

        !> Integer workspace
        integer(FP_SIZE), allocatable :: iwrk(:)

    contains

        !> MSE accessor (shared by all types)
        procedure, non_overridable :: mse => fitter_mse

        !> Base field helpers for comm (non-overridable, called by subtypes)
        procedure, non_overridable :: core_comm_size   => fitter_core_comm_size
        procedure, non_overridable :: core_comm_pack   => fitter_core_comm_pack
        procedure, non_overridable :: core_comm_expand => fitter_core_comm_expand

        !> Base field reset (called by subtype destroy methods)
        procedure, non_overridable :: destroy_base => fitter_destroy_base

        !> Deferred communication interface
        procedure(comm_size_if),   deferred :: comm_size
        procedure(comm_pack_if),   deferred :: comm_pack
        procedure(comm_expand_if), deferred :: comm_expand

    end type fitpack_fitter

    abstract interface
        elemental integer(FP_SIZE) function comm_size_if(this)
            import fitpack_fitter, FP_SIZE
            class(fitpack_fitter), intent(in) :: this
        end function
        pure subroutine comm_pack_if(this, buffer)
            import fitpack_fitter, FP_COMM
            class(fitpack_fitter), intent(in) :: this
            real(FP_COMM), intent(out) :: buffer(:)
        end subroutine
        pure subroutine comm_expand_if(this, buffer)
            import fitpack_fitter, FP_COMM
            class(fitpack_fitter), intent(inout) :: this
            real(FP_COMM), intent(in) :: buffer(:)
        end subroutine
    end interface

contains

    !> Return fitting MSE
    elemental real(FP_REAL) function fitter_mse(this)
        class(fitpack_fitter), intent(in) :: this
        fitter_mse = this%fp
    end function fitter_mse

    !> Number of FP_COMM slots needed for base fields:
    !> iopt (1 scalar) + smoothing (1 scalar) + fp (1 scalar) + c(:) + iwrk(:)
    elemental integer(FP_SIZE) function fitter_core_comm_size(this)
        class(fitpack_fitter), intent(in) :: this
        fitter_core_comm_size = 3 + FP_COMM_SIZE(this%c) + FP_COMM_SIZE(this%iwrk)
    end function fitter_core_comm_size

    !> Pack base fields into communication buffer
    pure subroutine fitter_core_comm_pack(this, buffer)
        class(fitpack_fitter), intent(in) :: this
        real(FP_COMM), intent(out) :: buffer(:)

        integer(FP_SIZE) :: pos

        buffer(1) = real(this%iopt, FP_COMM)
        buffer(2) = this%smoothing
        buffer(3) = this%fp
        pos = 4
        call FP_COMM_PACK(this%c, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%c)
        call FP_COMM_PACK(this%iwrk, buffer(pos:))

    end subroutine fitter_core_comm_pack

    !> Expand base fields from communication buffer
    pure subroutine fitter_core_comm_expand(this, buffer)
        class(fitpack_fitter), intent(inout) :: this
        real(FP_COMM), intent(in) :: buffer(:)

        integer(FP_SIZE) :: pos

        this%iopt      = nint(buffer(1), FP_FLAG)
        this%smoothing = buffer(2)
        this%fp        = buffer(3)
        pos = 4
        call FP_COMM_EXPAND(this%c, buffer(pos:));    pos = pos + FP_COMM_SIZE(this%c)
        call FP_COMM_EXPAND(this%iwrk, buffer(pos:))

    end subroutine fitter_core_comm_expand

    !> Reset base fields to defaults and deallocate arrays
    elemental subroutine fitter_destroy_base(this)
        class(fitpack_fitter), intent(inout) :: this
        integer :: ierr
        this%iopt      = IOPT_NEW_SMOOTHING
        this%smoothing = 1000.0_FP_REAL
        this%fp        = zero
        deallocate(this%c, stat=ierr)
        deallocate(this%iwrk, stat=ierr)
    end subroutine fitter_destroy_base

end module fitpack_fitters
