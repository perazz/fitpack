    // ===========================================================================================
    // extra_methods/fpCurveInherit.hpp (hand-maintained)
    //
    // This class overrides the raw new_fit / fit / interpolate / eval, which hides the
    // ergonomic overloads spliced into fpCurve. Re-expose them; the overrides above still win
    // for their own signatures, and every ergonomic overload dispatches back through them.
    // ===========================================================================================

    using fpCurve::new_fit;
    using fpCurve::fit;
    using fpCurve::interpolate;
    using fpCurve::eval;
