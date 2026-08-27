    // ===========================================================================================
    // fpPoint<dim> overloads — extra_methods/fpGridSplinePoints.hpp (hand-maintained)
    //
    // The N-D analog of the parametric-curve point overloads: a scattered evaluation site of an
    // N-dimensional gridded spline is a point, and a std::vector<fpPoint<dim>> is bit-identical
    // to the Fortran xp(dim,np) it feeds. `dim` is deduced from the argument; it is checked
    // against the spline's own dims() and reports FITPACK_INPUT_ERROR on a mismatch.
    // ===========================================================================================

    //! @brief Evaluate the spline at one scattered point.
    template <FP_SIZE dim>
    FP_REAL eval(const fpPoint<dim>& x, FP_FLAG* ierr = nullptr) const
    {
        if (dims() != dim) { if (ierr) *ierr = FITPACK_INPUT_ERROR; return FP_REAL(0); }
        std::vector<FP_REAL> xw(x.begin(), x.end());
        return eval(xw, ierr);
    }

    //! @brief Evaluate the spline at every scattered point.
    template <FP_SIZE dim>
    std::vector<FP_REAL> eval(const std::vector<fpPoint<dim>>& xp, FP_FLAG* ierr = nullptr)
    {
        if (dims() != dim) { if (ierr) *ierr = FITPACK_INPUT_ERROR; return {}; }
        std::vector<FP_REAL> xw = fpPointFlatten<dim>(xp);
        return eval(dim, static_cast<FP_SIZE>(xp.size()), xw, ierr);
    }

    //! @brief Partial derivative of orders nu(1..dim) at one scattered point.
    template <FP_SIZE dim>
    FP_REAL dfdx(const fpPoint<dim>& x, const std::vector<FP_SIZE>& nu, FP_FLAG* ierr = nullptr) const
    {
        if (dims() != dim) { if (ierr) *ierr = FITPACK_INPUT_ERROR; return FP_REAL(0); }
        std::vector<FP_REAL> xw(x.begin(), x.end());
        std::vector<FP_SIZE> nuw(nu);
        return dfdx(xw, nuw, ierr);
    }

    //! @brief Partial derivative of orders nu(1..dim) at every scattered point.
    template <FP_SIZE dim>
    std::vector<FP_REAL> dfdx(const std::vector<fpPoint<dim>>& xp, const std::vector<FP_SIZE>& nu,
                              FP_FLAG* ierr = nullptr)
    {
        if (dims() != dim) { if (ierr) *ierr = FITPACK_INPUT_ERROR; return {}; }
        std::vector<FP_REAL> xw = fpPointFlatten<dim>(xp);
        std::vector<FP_SIZE> nuw(nu);
        return dfdx(dim, static_cast<FP_SIZE>(xp.size()), xw, nuw, ierr);
    }

    //! @brief Integral of the spline over the box [lower, upper].
    template <FP_SIZE dim>
    FP_REAL integral(const fpPoint<dim>& lower, const fpPoint<dim>& upper) const
    {
        if (dims() != dim) return FP_REAL(0);
        std::vector<FP_REAL> lw(lower.begin(), lower.end()), uw(upper.begin(), upper.end());
        return integral(lw, uw);
    }
