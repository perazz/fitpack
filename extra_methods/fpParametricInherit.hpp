    // ===========================================================================================
    // extra_methods/fpParametricInherit.hpp (hand-maintained)
    //
    // This class overrides the raw new_fit / fit, which hides the fpPoint<dim> overloads
    // spliced into fpParametricCurve. Re-expose them; the overrides above still win for their
    // own signatures, and every fpPoint<dim> overload dispatches back through them.
    // ===========================================================================================

    using fpParametricCurve::new_fit;
    using fpParametricCurve::fit;
