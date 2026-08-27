    // ===========================================================================================
    // fpPoint<dim> overloads — extra_methods/fpParametricPoints.hpp (hand-maintained)
    //
    // The dimension is a template argument, so a point is a fixed-size POD rather than a heap
    // allocation, and a std::vector<fpPoint<dim>> is bit-identical to the Fortran x(dim,m) it
    // feeds. Generalizes the idea of PR #61 (@illionj) to the whole parametric family.
    //
    // Everything here forwards to a virtual raw method, so fpClosedCurve and fpConstrainedCurve
    // inherit these overloads and still reach their own Fortran routines through dispatch.
    // Methods that read an existing fit check idim() against dim and report
    // FITPACK_INPUT_ERROR on a mismatch rather than reading past the result buffer.
    // ===========================================================================================

    //! @brief Spline degree k of the current fit.
    FP_SIZE degree() const { return order(); }

    //! @brief Number of space dimensions of the fitted curve.
    FP_SIZE ndim() const { return idim(); }

    //! @brief Fit a new curve through the points, choosing u by cumulative Euclidean distance.
    template <FP_SIZE dim>
    FP_FLAG new_fit(const std::vector<fpPoint<dim>>& x, FP_REAL smoothing = 1000.0, FP_SIZE order = 3)
    {
        std::vector<FP_REAL> xw = fpPointFlatten<dim>(x);
        return new_fit(dim, static_cast<FP_SIZE>(x.size()), xw, nullptr, nullptr, &smoothing, &order);
    }

    //! @brief Fit a new curve through the points at the given parameter values u.
    template <FP_SIZE dim>
    FP_FLAG new_fit(const std::vector<fpPoint<dim>>& x, const std::vector<FP_REAL>& u,
                    FP_REAL smoothing = 1000.0, FP_SIZE order = 3)
    {
        std::vector<FP_REAL> xw = fpPointFlatten<dim>(x), uw(u);
        return new_fit(dim, static_cast<FP_SIZE>(x.size()), xw, uw.data(), nullptr, &smoothing, &order);
    }

    //! @brief Fit a new curve through the points at parameter values u, with weights w.
    template <FP_SIZE dim>
    FP_FLAG new_fit(const std::vector<fpPoint<dim>>& x, const std::vector<FP_REAL>& u,
                    const std::vector<FP_REAL>& w, FP_REAL smoothing = 1000.0, FP_SIZE order = 3)
    {
        std::vector<FP_REAL> xw = fpPointFlatten<dim>(x), uw(u), ww(w);
        return new_fit(dim, static_cast<FP_SIZE>(x.size()), xw, uw.data(), ww.data(), &smoothing, &order);
    }

    //! @brief Refit the current points with a new spline order.
    FP_FLAG fit(FP_SIZE order)                    { return fit(nullptr, &order, nullptr); }
    //! @brief Refit the current points with a new smoothing.
    FP_FLAG fit(FP_REAL smoothing)                { return fit(&smoothing, nullptr, nullptr); }
    //! @brief Refit the current points with a new smoothing and spline order.
    FP_FLAG fit(FP_REAL smoothing, FP_SIZE order) { return fit(&smoothing, &order, nullptr); }

    //! @brief Evaluate the curve at the parameter value u.
    template <FP_SIZE dim>
    fpPoint<dim> eval(FP_REAL u, FP_FLAG* ierr = nullptr)
    {
        fpPoint<dim> y{};
        if (idim() != dim) { if (ierr) *ierr = FITPACK_INPUT_ERROR; return y; }

        FP_FLAG ierr0 = FITPACK_OK;
        const std::vector<FP_REAL> flat = eval_one(u, dim, &ierr0);
        if (ierr) *ierr = ierr0;
        for (FP_SIZE j = 0; j < dim; ++j) y[j] = flat[static_cast<std::size_t>(j)];
        return y;
    }

    //! @brief Evaluate the curve at every parameter value in u.
    template <FP_SIZE dim>
    std::vector<fpPoint<dim>> eval(const std::vector<FP_REAL>& u, FP_FLAG* ierr = nullptr)
    {
        if (idim() != dim) { if (ierr) *ierr = FITPACK_INPUT_ERROR; return {}; }

        std::vector<FP_REAL> uw(u);
        FP_FLAG ierr0 = FITPACK_OK;
        const std::vector<FP_REAL> flat = eval_many(uw, dim * static_cast<FP_SIZE>(u.size()), &ierr0);
        if (ierr) *ierr = ierr0;
        return fpPointGather<dim>(flat);
    }

    //! @brief Derivative of the given order at the parameter value u.
    template <FP_SIZE dim>
    fpPoint<dim> ddu(FP_REAL u, FP_SIZE order, FP_FLAG* ierr = nullptr)
    {
        fpPoint<dim> d{};
        if (idim() != dim) { if (ierr) *ierr = FITPACK_INPUT_ERROR; return d; }

        FP_FLAG ierr0 = FITPACK_OK;
        const std::vector<FP_REAL> flat = dfdx(u, order, dim, &ierr0);
        if (ierr) *ierr = ierr0;
        for (FP_SIZE j = 0; j < dim; ++j) d[j] = flat[static_cast<std::size_t>(j)];
        return d;
    }

    //! @brief All derivatives (orders 0..k) at the parameter value u.
    template <FP_SIZE dim>
    std::vector<fpPoint<dim>> ddu_all(FP_REAL u, FP_FLAG* ierr = nullptr)
    {
        if (idim() != dim) { if (ierr) *ierr = FITPACK_INPUT_ERROR; return {}; }

        FP_FLAG ierr0 = FITPACK_OK;
        const std::vector<FP_REAL> flat = dfdx_all(u, dim * (degree() + 1), &ierr0);
        if (ierr) *ierr = ierr0;
        return fpPointGather<dim>(flat);
    }
