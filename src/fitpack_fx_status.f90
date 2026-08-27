!   ***********************************************************************************************
!   **                                         FITPACK                                         **
!   **                     Modern Fortran Fitting Package — C/C++ Bindings                     **
!   ***********************************************************************************************
!   **    fitpack_fx_status
!> @brief Self-contained fx_status shim for standalone bindings (no fortran-arrays dependency)
!   ***********************************************************************************************
!> @author Binding Generator
!   ***********************************************************************************************
!
!   Layout contract. type(fx_status) below is the Fortran half of the C struct that the
!   generated `*_c.h` headers inline (templates/c_header_standalone.h.jinja2):
!
!       #define FX_LEN_STATUS_MSG 248
!       typedef struct fx_status { bool ok; int code; char message[FX_LEN_STATUS_MSG]; } fx_status;
!
!   Field order, kinds and LEN_STATUS_MSG == FX_LEN_STATUS_MSG must stay in lock-step with
!   that struct AND with fortran-arrays' own arrays_c definition, so a consumer can link
!   standalone wrappers and real arrays_c wrappers into the same binary.
!
module fitpack_fx_status
    use, intrinsic :: iso_c_binding
    implicit none(type, external)
    private

    public :: fx_status
    public :: fx_status_ok
    public :: handle_error
    public :: FX_OK, FX_ERROR_ALLOCATION, FX_ERROR_DEALLOCATION

    integer, parameter :: LEN_STATUS_MSG = 248

    !> Status codes. Values mirror arrays_c; only the codes the generated
    !> wrappers actually raise are carried here.
    integer(c_int), parameter :: FX_OK = 0
    integer(c_int), parameter :: FX_ERROR_ALLOCATION = 1
    integer(c_int), parameter :: FX_ERROR_DEALLOCATION = 2

    !> C-interoperable status type for error reporting across language boundaries.
    type, bind(C) :: fx_status
        logical(c_bool)   :: ok = .true._c_bool          ! Quick success check
        integer(c_int)    :: code = FX_OK                ! Error code for programmatic handling
        character(c_char) :: message(LEN_STATUS_MSG)     ! Human-readable error message
    end type fx_status

    !> Success status constant (initialization without a procedure call)
    type(fx_status), parameter :: fx_status_ok = fx_status( &
        ok = .true._c_bool, &
        code = FX_OK, &
        message = c_null_char)

    !> Error-status construction: stat0 = fx_status(FX_ERROR_ALLOCATION, "allocation failed")
    interface fx_status
        module procedure fx_status_new
    end interface fx_status

    contains

    !> Build an error status from a code and a Fortran string
    pure function fx_status_new(code, message) result(stat)
        integer(c_int), intent(in) :: code
        character(*), intent(in) :: message
        type(fx_status) :: stat

        integer :: i, n

        stat%ok = .false._c_bool
        stat%code = code

        n = min(len_trim(message), LEN_STATUS_MSG - 1)
        do i = 1, n
            stat%message(i) = message(i:i)
        end do
        stat%message(n + 1) = c_null_char
    end function fx_status_new

    !> Read a status message back as a Fortran string
    pure function fx_status_message(status) result(msg)
        type(fx_status), intent(in) :: status
        character(len=LEN_STATUS_MSG) :: msg
        integer :: i

        msg = ''
        do i = 1, LEN_STATUS_MSG
            if (status%message(i) == c_null_char) exit
            msg(i:i) = status%message(i)
        end do
    end function fx_status_message

    !> Copy the local status to the caller's optional out-argument, or stop on error
    !> @param stat0 - Local status carrying the error information
    !> @param stat  - Optional output status for the caller
    subroutine handle_error(stat0, stat)
        type(fx_status), intent(in) :: stat0
        type(fx_status), intent(out), optional :: stat

        if (present(stat)) then
            stat = stat0
        elseif (.not.stat0%ok) then
            error stop fx_status_message(stat0)
        end if
    end subroutine handle_error

end module fitpack_fx_status
