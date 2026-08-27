    // ===========================================================================================
    // Ergonomic overloads — extra_methods/fpCurveErgonomics.hpp (hand-maintained)
    //
    // The call shapes of the pre-0.3.0 hand-written wrapper, on top of the generated raw entry
    // points. Every one forwards to a virtual raw method, so fpPeriodicCurve inherits them and
    // still reaches its own Fortran routines through dispatch.
    // ===========================================================================================

    //! @brief Spline degree k of the current fit.
    FP_SIZE degree() const { return order(); }

    //! @brief Fit a new curve through (x, y).
    FP_FLAG new_fit(const std::vector<FP_REAL>& x, const std::vector<FP_REAL>& y,
                    FP_REAL smoothing = 1000.0)
    {
        std::vector<FP_REAL> xw(x), yw(y);
        return new_fit(xw, yw, nullptr, &smoothing, nullptr);
    }

    //! @brief Fit a new curve through (x, y) with per-point weights w.
    FP_FLAG new_fit(const std::vector<FP_REAL>& x, const std::vector<FP_REAL>& y,
                    const std::vector<FP_REAL>& w, FP_REAL smoothing = 1000.0)
    {
        std::vector<FP_REAL> xw(x), yw(y), ww(w);
        return new_fit(xw, yw, ww.data(), &smoothing, nullptr);
    }

    //! @brief Refit the current points with a new spline order.
    FP_FLAG fit(FP_SIZE order)                    { return fit(nullptr, &order, nullptr); }
    //! @brief Refit the current points with a new smoothing.
    FP_FLAG fit(FP_REAL smoothing)                { return fit(&smoothing, nullptr, nullptr); }
    //! @brief Refit the current points with a new smoothing and spline order.
    FP_FLAG fit(FP_REAL smoothing, FP_SIZE order) { return fit(&smoothing, &order, nullptr); }

    //! @brief Interpolating fit of the given spline order.
    FP_FLAG interpolate(FP_SIZE order) { return interpolate(&order, nullptr); }

    //! @brief Evaluate the curve at x, reporting the status through @p ierr.
    FP_REAL eval(FP_REAL x, FP_FLAG* ierr)
    {
        FP_FLAG ierr0 = FITPACK_OK;
        const FP_REAL y = eval(x, ierr0);
        if (ierr) *ierr = ierr0;
        return y;
    }

    //! @brief Evaluate the curve at every x, reporting the status through @p ierr.
    std::vector<FP_REAL> eval(const std::vector<FP_REAL>& x, FP_FLAG* ierr)
    {
        std::vector<FP_REAL> xw(x);
        FP_FLAG ierr0 = FITPACK_OK;
        std::vector<FP_REAL> y = eval(xw, ierr0);
        if (ierr) *ierr = ierr0;
        return y;
    }

    //! @brief Evaluate the curve at @p npts points evenly spaced over [xmin, xmax].
    //! @return The (x, y) pairs, so the result plots without a second array.
    //! @note @p npts has no default: `eval(x, ierr)` with an integer `ierr` would otherwise be
    //!       ambiguous against the scalar `eval` under ISO C++ (gcc accepts it, clang does not).
    std::vector<fpPoint<2>> eval(FP_REAL xmin, FP_REAL xmax, FP_SIZE npts, FP_FLAG* ierr = nullptr)
    {
        std::vector<FP_REAL> x(static_cast<std::size_t>(npts > 0 ? npts : 0));
        if (npts > 1)
        {
            const FP_REAL dx = (xmax - xmin) / (npts - 1);
            for (FP_SIZE i = 0; i < npts; ++i) x[static_cast<std::size_t>(i)] = xmin + i * dx;
        }
        else if (npts == 1) x[0] = xmin;

        FP_FLAG ierr0 = FITPACK_OK;
        std::vector<FP_REAL> y = eval(x, ierr0);
        if (ierr) *ierr = ierr0;

        std::vector<fpPoint<2>> xy(x.size());
        for (std::size_t i = 0; i < x.size() && i < y.size(); ++i) xy[i] = {x[i], y[i]};
        return xy;
    }

    //! @brief Derivative of the given order at x.
    //! @note Not const: the Fortran receiver of curve_derivative is intent(inout).
    FP_REAL ddx(FP_REAL x, FP_SIZE order, FP_FLAG* ierr = nullptr) { return dfdx(x, order, ierr); }

    //! @brief All derivatives (orders 0..k) at x.
    std::vector<FP_REAL> ddx(FP_REAL x, FP_FLAG* ierr = nullptr)
    {
        FP_FLAG ierr0 = FITPACK_OK;
        std::vector<FP_REAL> d = dfdx_all(x, ierr0, degree() + 1);
        if (ierr) *ierr = ierr0;
        return d;
    }

    //! @brief Spline behaviour outside the support (one of the OUTSIDE_* flags).
    void    set_bc(FP_FLAG value) { bc() = value; }
    FP_FLAG get_bc() const        { return bc(); }

    //! @brief Fourier coefficients A, B of the curve at the given angular frequencies.
    //! @return The FITPACK status flag; A and B are resized to match @p alpha.
    FP_FLAG fourier(const std::vector<FP_REAL>& alpha, std::vector<FP_REAL>& A, std::vector<FP_REAL>& B)
    {
        std::vector<FP_REAL> aw(alpha);
        A.resize(alpha.size());
        B.resize(alpha.size());
        FP_FLAG ierr = FITPACK_OK;
        fourier_coefficients(aw, A, B, &ierr);
        return ierr;
    }
