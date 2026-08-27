    // ===========================================================================================
    // Endpoint constraints — extra_methods/fpConstrainedPoints.hpp (hand-maintained)
    //
    // Each list holds the position and the successive derivatives to pin at that endpoint:
    // entry j is the j-th derivative, so a single entry fixes the point only.
    //
    // The const_cast is safe: the Fortran dummy is `optional, intent(in)`. The generated
    // wrapper takes a non-const pointer only because an OPTIONAL argument is passed by
    // address, and nullptr is how a C caller says "absent".
    // ===========================================================================================

    //! @brief Pin the begin endpoint; leave the end endpoint free.
    template <FP_SIZE dim>
    FP_FLAG constrain_begin(const std::vector<fpPoint<dim>>& ddx_begin)
    {
        FP_FLAG ierr = FITPACK_OK;
        set_constraints(dim, static_cast<FP_SIZE>(ddx_begin.size()), dim, 0,
                        const_cast<FP_REAL*>(fpPointData<dim>(ddx_begin)), nullptr, &ierr);
        return ierr;
    }

    //! @brief Pin both endpoints.
    template <FP_SIZE dim>
    FP_FLAG constrain_both(const std::vector<fpPoint<dim>>& ddx_begin,
                           const std::vector<fpPoint<dim>>& ddx_end)
    {
        FP_FLAG ierr = FITPACK_OK;
        set_constraints(dim, static_cast<FP_SIZE>(ddx_begin.size()),
                        dim, static_cast<FP_SIZE>(ddx_end.size()),
                        const_cast<FP_REAL*>(fpPointData<dim>(ddx_begin)),
                        const_cast<FP_REAL*>(fpPointData<dim>(ddx_end)), &ierr);
        return ierr;
    }

    //! @brief Pin the end endpoint; leave the begin endpoint free.
    template <FP_SIZE dim>
    FP_FLAG constrain_end(const std::vector<fpPoint<dim>>& ddx_end)
    {
        FP_FLAG ierr = FITPACK_OK;
        set_constraints(dim, 0, dim, static_cast<FP_SIZE>(ddx_end.size()),
                        nullptr, const_cast<FP_REAL*>(fpPointData<dim>(ddx_end)), &ierr);
        return ierr;
    }
