/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fpConstrainedCurve.hpp (class fpConstrainedCurve)
!> @brief Standalone C++ wrapper for fitpack_constrained_curve (no fortran-arrays dependency)
!
!   @author Federico Perini
!   @date   2026-08-27
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************/

#ifndef FPCONSTRAINEDCURVE_HPP_INCLUDED
#define FPCONSTRAINEDCURVE_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_constrained_curve_c.h"
#include "fpParametricCurve.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <variant>
#include <optional>
#include <array>

static_assert(sizeof(fitpack_constrained_curve_c) == sizeof(fitpack_parametric_curve_c),
    "C descriptor layout mismatch: fitpack_constrained_curve_c vs fitpack_parametric_curve_c");

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_constrained_curve
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_constrained_curve, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 * Extends fpParametricCurve (Fortran: extends(fitpack_parametric_curve))
 */
class fpConstrainedCurve : public fpParametricCurve {
public:
    // ===========================================================================================
    // Constructors and Destructor
    // ===========================================================================================

    /**
     * @brief Default constructor - allocates new fitpack_constrained_curve
     */
    fpConstrainedCurve() : fpParametricCurve(NoAlloc{}) {
        fitpack_constrained_curve_c_allocate(as<fitpack_constrained_curve_c>(), nullptr);
    }

    /**
     * @brief Destructor - deallocates if owned
     */
    ~fpConstrainedCurve() override {
        fitpack_constrained_curve_c_destroy(as<fitpack_constrained_curve_c>(), nullptr);
    }

    /**
     * @brief Copy constructor - deep copy
     */
    fpConstrainedCurve(const fpConstrainedCurve& other) : fpParametricCurve(NoAlloc{}) {
        *as<fitpack_constrained_curve_c>() = fitpack_constrained_curve_c_null;
        fitpack_constrained_curve_c_copy(as<fitpack_constrained_curve_c>(), other.as<fitpack_constrained_curve_c>(), false, nullptr);
    }

    /**
     * @brief Copy assignment - deep copy
     */
    fpConstrainedCurve& operator=(const fpConstrainedCurve& other) {
        if (this != &other) {
            fitpack_constrained_curve_c_destroy(as<fitpack_constrained_curve_c>(), nullptr);
            fitpack_constrained_curve_c_copy(as<fitpack_constrained_curve_c>(), other.as<fitpack_constrained_curve_c>(), false, nullptr);
        }
        return *this;
    }

    /**
     * @brief Move constructor - transfer ownership
     */
    fpConstrainedCurve(fpConstrainedCurve&& other) noexcept : fpParametricCurve(NoAlloc{}) {
        *as<fitpack_constrained_curve_c>() = fitpack_constrained_curve_c_null;
        fitpack_constrained_curve_c_move_alloc(as<fitpack_constrained_curve_c>(), other.as<fitpack_constrained_curve_c>(), nullptr);
    }

    /**
     * @brief Move assignment - transfer ownership
     */
    fpConstrainedCurve& operator=(fpConstrainedCurve&& other) noexcept {
        if (this != &other) {
            fitpack_constrained_curve_c_destroy(as<fitpack_constrained_curve_c>(), nullptr);
            fitpack_constrained_curve_c_move_alloc(as<fitpack_constrained_curve_c>(), other.as<fitpack_constrained_curve_c>(), nullptr);
        }
        return *this;
    }

    /**
     * @brief Construct from existing C wrapper (takes ownership if move=true)
     */
    explicit fpConstrainedCurve(fitpack_constrained_curve_c& c_wrapper, bool move = false) : fpParametricCurve(NoAlloc{}) {
        if (move) {
            fitpack_constrained_curve_c_move_alloc(as<fitpack_constrained_curve_c>(), &c_wrapper, nullptr);
        } else {
            *as<fitpack_constrained_curve_c>() = fitpack_constrained_curve_c_null;
            fitpack_constrained_curve_c_copy(as<fitpack_constrained_curve_c>(), &c_wrapper, false, nullptr);
        }
    }

    /**
     * @brief Non-owning view ctor: bit-copy the C handle, mark
     * is_pointer=true. Used by parent classes' polymorphic accessors
     * to wrap a base-class handle as the correct concrete subtype.
     */
    explicit fpConstrainedCurve(fitpack_constrained_curve_c& c_wrapper, ViewTag) : fpParametricCurve(NoAlloc{}) {
        *as<fitpack_constrained_curve_c>() = c_wrapper;
        as<fitpack_constrained_curve_c>()->is_pointer = true;
    }

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    bool is_allocated() const override {
        return as<fitpack_constrained_curve_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    bool is_pointer() const override {
        return as<fitpack_constrained_curve_c>()->is_pointer;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_constrained_curve_c").
     */
    const char* c_type_name() const override {
        return fitpack_constrained_curve_c_c_type_name(*as<fitpack_constrained_curve_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    const char* cpp_type_name() const override {
        return "fpConstrainedCurve";
    }

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_constrained_curve_c& c_handle() {
        return *as<fitpack_constrained_curve_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_constrained_curve_c& c_handle() const {
        return *as<fitpack_constrained_curve_c>();
    }

    /**
     * @brief Construct an empty (non-owning) wrapper, for use as a view target.
     *
     * The returned object has a null handle (cptr==nullptr, is_pointer==false)
     * and does NOT allocate a Fortran object. Intended as a fill target for
     * member-view accessors emitted on parent types — the parent populates
     * the handle with c_loc() into its own storage and sets is_pointer=true.
     */
    static fpConstrainedCurve make_view() {
        return fpConstrainedCurve(NoAlloc{});
    }

    /**
     * @brief Upcast to parent type (reference, no copy)
     */
    fpParametricCurve& as_parent() {
        return static_cast<fpParametricCurve&>(*this);
    }
    const fpParametricCurve& as_parent() const {
        return static_cast<const fpParametricCurve&>(*this);
    }

    // ===========================================================================================
    // Method Wrappers (standalone — no fxArray dependency)
    // ===========================================================================================

    virtual void clean_constraints() {
        fitpack_constrained_curve_c_clean_constraints(as<fitpack_constrained_curve_c>());
    }

    /**
     * @brief set_constraints
     */
    virtual void set_constraints(int32_t ddx_begin_n1, int32_t ddx_begin_n2, int32_t ddx_end_n1, int32_t ddx_end_n2, double* ddx_begin = nullptr, double* ddx_end = nullptr, int32_t* ierr = nullptr) {
        fitpack_constrained_curve_c_set_constraints(as<fitpack_constrained_curve_c>(), ddx_begin_n1, ddx_begin_n2, ddx_begin, ddx_end_n1, ddx_end_n2, ddx_end, ierr);
    }


    int32_t comm_size() const override {
        return fitpack_constrained_curve_c_comm_size(as<fitpack_constrained_curve_c>());
    }

    /**
     * @brief comm_pack
     */
    void comm_pack(std::vector<double>& buffer) const override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_constrained_curve_c_comm_pack(as<fitpack_constrained_curve_c>(), n, buffer.data());
    }


    /**
     * @brief comm_expand
     */
    void comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_constrained_curve_c_comm_expand(as<fitpack_constrained_curve_c>(), n, buffer.data());
    }


    /**
     * @brief new_points
     */
    void new_points(int32_t x_n1, int32_t x_n2, std::vector<double>& x, double* u = nullptr, double* w = nullptr) override {
        fitpack_constrained_curve_c_new_points(as<fitpack_constrained_curve_c>(), x_n1, x_n2, x.data(), u, w);
    }


    void set_default_parameters() override {
        fitpack_constrained_curve_c_set_default_parameters(as<fitpack_constrained_curve_c>());
    }

    /**
     * @brief new_fit
     */
    int32_t new_fit(int32_t x_n1, int32_t x_n2, std::vector<double>& x, double* u = nullptr, double* w = nullptr, double* smoothing = nullptr, int32_t* order = nullptr) override {
        return fitpack_constrained_curve_c_new_fit(as<fitpack_constrained_curve_c>(), x_n1, x_n2, x.data(), u, w, smoothing, order);
    }


    int32_t fit(double* smoothing = nullptr, int32_t* order = nullptr, bool* keep_knots = nullptr) override {
        return fitpack_constrained_curve_c_fit(as<fitpack_constrained_curve_c>(), smoothing, order, keep_knots);
    }

    int32_t interpolate(int32_t* order = nullptr, bool* reset_knots = nullptr) override {
        return fitpack_constrained_curve_c_interpolate(as<fitpack_constrained_curve_c>(), order, reset_knots);
    }

    int32_t least_squares(double* smoothing = nullptr, bool* reset_knots = nullptr) override {
        return fitpack_constrained_curve_c_least_squares(as<fitpack_constrained_curve_c>(), smoothing, reset_knots);
    }

    /**
     * @brief eval_one
     */
    std::vector<double> eval_one(double u, int32_t n_result, int32_t* ierr = nullptr) override {
        std::vector<double> result(n_result);
        fitpack_constrained_curve_c_eval_one(as<fitpack_constrained_curve_c>(), u, ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief eval_many
     */
    std::vector<double> eval_many(std::vector<double>& u, int32_t n_result, int32_t* ierr = nullptr) override {
        int32_t n = static_cast<int32_t>(u.size());
        std::vector<double> result(n_result);
        fitpack_constrained_curve_c_eval_many(as<fitpack_constrained_curve_c>(), n, u.data(), ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx
     */
    std::vector<double> dfdx(double u, int32_t order, int32_t n_result, int32_t* ierr = nullptr) override {
        std::vector<double> result(n_result);
        fitpack_constrained_curve_c_curve_derivative(as<fitpack_constrained_curve_c>(), u, order, ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx
     */
    std::vector<double> dfdx(std::vector<double>& u, int32_t order, int32_t n_result, int32_t* ierr = nullptr) override {
        int32_t n = static_cast<int32_t>(u.size());
        std::vector<double> result(n_result);
        fitpack_constrained_curve_c_curve_derivatives(as<fitpack_constrained_curve_c>(), n, u.data(), order, ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx_all
     */
    std::vector<double> dfdx_all(double u, int32_t n_result, int32_t* ierr = nullptr) override {
        std::vector<double> result(n_result);
        fitpack_constrained_curve_c_curve_all_derivatives(as<fitpack_constrained_curve_c>(), u, ierr, result.data(), n_result);
        return result;
    }


    double mse() const override {
        return fitpack_constrained_curve_c_mse(as<fitpack_constrained_curve_c>());
    }

    int32_t core_comm_size() const override {
        return fitpack_constrained_curve_c_core_comm_size(as<fitpack_constrained_curve_c>());
    }

    /**
     * @brief core_comm_pack
     */
    void core_comm_pack(std::vector<double>& buffer) const override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_constrained_curve_c_core_comm_pack(as<fitpack_constrained_curve_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    void core_comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_constrained_curve_c_core_comm_expand(as<fitpack_constrained_curve_c>(), n, buffer.data());
    }


    void destroy_base() override {
        fitpack_constrained_curve_c_destroy_base(as<fitpack_constrained_curve_c>());
    }

    // ===========================================================================================
    // Component Array Accessors
    // ===========================================================================================

    /**
     * @brief Deep copy of component 'deriv_begin' as a std::vector.
     *
     * The rank-2 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with deriv_begin_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> deriv_begin_vector() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_deriv_begin(as<fitpack_constrained_curve_c>(), &raw, extents);
        const int64_t total = extents[0] * extents[1];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

    /**
     * @brief Extents of component 'deriv_begin', leading dimension first.
     * @return All zeros when the component is unallocated.
     */
    std::array<int64_t, 2> deriv_begin_shape() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_deriv_begin(as<fitpack_constrained_curve_c>(), &raw, extents);
        return { extents[0], extents[1] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'deriv_begin'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through deriv_begin()(i, j) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> deriv_begin() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_deriv_begin(as<fitpack_constrained_curve_c>(), &raw, extents);
        FX_SIZE bounds[4];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 2; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "deriv_begin", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(2), bounds);
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    /**
     * @brief Deep copy of component 'deriv_end' as a std::vector.
     *
     * The rank-2 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with deriv_end_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> deriv_end_vector() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_deriv_end(as<fitpack_constrained_curve_c>(), &raw, extents);
        const int64_t total = extents[0] * extents[1];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

    /**
     * @brief Extents of component 'deriv_end', leading dimension first.
     * @return All zeros when the component is unallocated.
     */
    std::array<int64_t, 2> deriv_end_shape() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_deriv_end(as<fitpack_constrained_curve_c>(), &raw, extents);
        return { extents[0], extents[1] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'deriv_end'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through deriv_end()(i, j) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> deriv_end() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_deriv_end(as<fitpack_constrained_curve_c>(), &raw, extents);
        FX_SIZE bounds[4];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 2; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "deriv_end", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(2), bounds);
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    /**
     * @brief Deep copy of component 'xx' as a std::vector.
     *
     * The rank-2 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with xx_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> xx_vector() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_xx(as<fitpack_constrained_curve_c>(), &raw, extents);
        const int64_t total = extents[0] * extents[1];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

    /**
     * @brief Extents of component 'xx', leading dimension first.
     * @return All zeros when the component is unallocated.
     */
    std::array<int64_t, 2> xx_shape() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_xx(as<fitpack_constrained_curve_c>(), &raw, extents);
        return { extents[0], extents[1] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'xx'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through xx()(i, j) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> xx() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_xx(as<fitpack_constrained_curve_c>(), &raw, extents);
        FX_SIZE bounds[4];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 2; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "xx", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(2), bounds);
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    /**
     * @brief Deep copy of component 'cp' as a std::vector.
     *
     * The rank-2 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with cp_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> cp_vector() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_cp(as<fitpack_constrained_curve_c>(), &raw, extents);
        const int64_t total = extents[0] * extents[1];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

    /**
     * @brief Extents of component 'cp', leading dimension first.
     * @return All zeros when the component is unallocated.
     */
    std::array<int64_t, 2> cp_shape() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_cp(as<fitpack_constrained_curve_c>(), &raw, extents);
        return { extents[0], extents[1] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'cp'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through cp()(i, j) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> cp() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_constrained_curve_c_getcomp_cp(as<fitpack_constrained_curve_c>(), &raw, extents);
        FX_SIZE bounds[4];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 2; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "cp", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(2), bounds);
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    // ===========================================================================================
    // Scalar Property Accessors
    // ===========================================================================================

    int32_t& ib() {
        return *fitpack_constrained_curve_c_ref_ib(as<fitpack_constrained_curve_c>());
    }
    const int32_t& ib() const {
        return *fitpack_constrained_curve_c_ref_ib(as<fitpack_constrained_curve_c>());
    }

    int32_t& ie() {
        return *fitpack_constrained_curve_c_ref_ie(as<fitpack_constrained_curve_c>());
    }
    const int32_t& ie() const {
        return *fitpack_constrained_curve_c_ref_ie(as<fitpack_constrained_curve_c>());
    }

    // ===========================================================================================
    // extra_methods/fpParametricInherit.hpp (hand-maintained)
    //
    // This class overrides the raw new_fit / fit, which hides the fpPoint<dim> overloads
    // spliced into fpParametricCurve. Re-expose them; the overrides above still win for their
    // own signatures, and every fpPoint<dim> overload dispatches back through them.
    // ===========================================================================================

    using fpParametricCurve::new_fit;
    using fpParametricCurve::fit;

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

protected:
    explicit fpConstrainedCurve(NoAlloc tag) : fpParametricCurve(tag) {}

};

#endif /* FPCONSTRAINEDCURVE_HPP_INCLUDED */
