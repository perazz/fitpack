/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fpClosedCurve.hpp (class fpClosedCurve)
!> @brief Standalone C++ wrapper for fitpack_closed_curve (no fortran-arrays dependency)
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

#ifndef FPCLOSEDCURVE_HPP_INCLUDED
#define FPCLOSEDCURVE_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_closed_curve_c.h"
#include "fpParametricCurve.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <variant>
#include <optional>

static_assert(sizeof(fitpack_closed_curve_c) == sizeof(fitpack_parametric_curve_c),
    "C descriptor layout mismatch: fitpack_closed_curve_c vs fitpack_parametric_curve_c");

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_closed_curve
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_closed_curve, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 * Extends fpParametricCurve (Fortran: extends(fitpack_parametric_curve))
 */
class fpClosedCurve : public fpParametricCurve {
public:
    // ===========================================================================================
    // Constructors and Destructor
    // ===========================================================================================

    /**
     * @brief Default constructor - allocates new fitpack_closed_curve
     */
    fpClosedCurve() : fpParametricCurve(NoAlloc{}) {
        fitpack_closed_curve_c_allocate(as<fitpack_closed_curve_c>(), nullptr);
    }

    /**
     * @brief Destructor - deallocates if owned
     */
    ~fpClosedCurve() override {
        fitpack_closed_curve_c_destroy(as<fitpack_closed_curve_c>(), nullptr);
    }

    /**
     * @brief Copy constructor - deep copy
     */
    fpClosedCurve(const fpClosedCurve& other) : fpParametricCurve(NoAlloc{}) {
        *as<fitpack_closed_curve_c>() = fitpack_closed_curve_c_null;
        fitpack_closed_curve_c_copy(as<fitpack_closed_curve_c>(), other.as<fitpack_closed_curve_c>(), false, nullptr);
    }

    /**
     * @brief Copy assignment - deep copy
     */
    fpClosedCurve& operator=(const fpClosedCurve& other) {
        if (this != &other) {
            fitpack_closed_curve_c_destroy(as<fitpack_closed_curve_c>(), nullptr);
            fitpack_closed_curve_c_copy(as<fitpack_closed_curve_c>(), other.as<fitpack_closed_curve_c>(), false, nullptr);
        }
        return *this;
    }

    /**
     * @brief Move constructor - transfer ownership
     */
    fpClosedCurve(fpClosedCurve&& other) noexcept : fpParametricCurve(NoAlloc{}) {
        *as<fitpack_closed_curve_c>() = fitpack_closed_curve_c_null;
        fitpack_closed_curve_c_move_alloc(as<fitpack_closed_curve_c>(), other.as<fitpack_closed_curve_c>(), nullptr);
    }

    /**
     * @brief Move assignment - transfer ownership
     */
    fpClosedCurve& operator=(fpClosedCurve&& other) noexcept {
        if (this != &other) {
            fitpack_closed_curve_c_destroy(as<fitpack_closed_curve_c>(), nullptr);
            fitpack_closed_curve_c_move_alloc(as<fitpack_closed_curve_c>(), other.as<fitpack_closed_curve_c>(), nullptr);
        }
        return *this;
    }

    /**
     * @brief Construct from existing C wrapper (takes ownership if move=true)
     */
    explicit fpClosedCurve(fitpack_closed_curve_c& c_wrapper, bool move = false) : fpParametricCurve(NoAlloc{}) {
        if (move) {
            fitpack_closed_curve_c_move_alloc(as<fitpack_closed_curve_c>(), &c_wrapper, nullptr);
        } else {
            *as<fitpack_closed_curve_c>() = fitpack_closed_curve_c_null;
            fitpack_closed_curve_c_copy(as<fitpack_closed_curve_c>(), &c_wrapper, false, nullptr);
        }
    }

    /**
     * @brief Non-owning view ctor: bit-copy the C handle, mark
     * is_pointer=true. Used by parent classes' polymorphic accessors
     * to wrap a base-class handle as the correct concrete subtype.
     */
    explicit fpClosedCurve(fitpack_closed_curve_c& c_wrapper, ViewTag) : fpParametricCurve(NoAlloc{}) {
        *as<fitpack_closed_curve_c>() = c_wrapper;
        as<fitpack_closed_curve_c>()->is_pointer = true;
    }

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    bool is_allocated() const override {
        return as<fitpack_closed_curve_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    bool is_pointer() const override {
        return as<fitpack_closed_curve_c>()->is_pointer;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_closed_curve_c").
     */
    const char* c_type_name() const override {
        return fitpack_closed_curve_c_c_type_name(*as<fitpack_closed_curve_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    const char* cpp_type_name() const override {
        return "fpClosedCurve";
    }

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_closed_curve_c& c_handle() {
        return *as<fitpack_closed_curve_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_closed_curve_c& c_handle() const {
        return *as<fitpack_closed_curve_c>();
    }

    /**
     * @brief Construct an empty (non-owning) wrapper, for use as a view target.
     *
     * The returned object has a null handle (cptr==nullptr, is_pointer==false)
     * and does NOT allocate a Fortran object. Intended as a fill target for
     * member-view accessors emitted on parent types — the parent populates
     * the handle with c_loc() into its own storage and sets is_pointer=true.
     */
    static fpClosedCurve make_view() {
        return fpClosedCurve(NoAlloc{});
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

    /**
     * @brief new_points
     */
    void new_points(int32_t x_n1, int32_t x_n2, std::vector<double>& x, double* u = nullptr, double* w = nullptr) override {
        fitpack_closed_curve_c_new_points(as<fitpack_closed_curve_c>(), x_n1, x_n2, x.data(), u, w);
    }


    void set_default_parameters() override {
        fitpack_closed_curve_c_set_default_parameters(as<fitpack_closed_curve_c>());
    }

    /**
     * @brief new_fit
     */
    int32_t new_fit(int32_t x_n1, int32_t x_n2, std::vector<double>& x, double* u = nullptr, double* w = nullptr, double* smoothing = nullptr, int32_t* order = nullptr) override {
        return fitpack_closed_curve_c_new_fit(as<fitpack_closed_curve_c>(), x_n1, x_n2, x.data(), u, w, smoothing, order);
    }


    int32_t fit(double* smoothing = nullptr, int32_t* order = nullptr, bool* keep_knots = nullptr) override {
        return fitpack_closed_curve_c_fit(as<fitpack_closed_curve_c>(), smoothing, order, keep_knots);
    }

    int32_t interpolate(int32_t* order = nullptr, bool* reset_knots = nullptr) override {
        return fitpack_closed_curve_c_interpolate(as<fitpack_closed_curve_c>(), order, reset_knots);
    }

    int32_t least_squares(double* smoothing = nullptr, bool* reset_knots = nullptr) override {
        return fitpack_closed_curve_c_least_squares(as<fitpack_closed_curve_c>(), smoothing, reset_knots);
    }

    /**
     * @brief eval_one
     */
    std::vector<double> eval_one(double u, int32_t n_result, int32_t* ierr = nullptr) override {
        std::vector<double> result(n_result);
        fitpack_closed_curve_c_eval_one(as<fitpack_closed_curve_c>(), u, ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief eval_many
     */
    std::vector<double> eval_many(std::vector<double>& u, int32_t n_result, int32_t* ierr = nullptr) override {
        int32_t n = static_cast<int32_t>(u.size());
        std::vector<double> result(n_result);
        fitpack_closed_curve_c_eval_many(as<fitpack_closed_curve_c>(), n, u.data(), ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx
     */
    std::vector<double> dfdx(double u, int32_t order, int32_t n_result, int32_t* ierr = nullptr) override {
        std::vector<double> result(n_result);
        fitpack_closed_curve_c_curve_derivative(as<fitpack_closed_curve_c>(), u, order, ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx
     */
    std::vector<double> dfdx(std::vector<double>& u, int32_t order, int32_t n_result, int32_t* ierr = nullptr) override {
        int32_t n = static_cast<int32_t>(u.size());
        std::vector<double> result(n_result);
        fitpack_closed_curve_c_curve_derivatives(as<fitpack_closed_curve_c>(), n, u.data(), order, ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx_all
     */
    std::vector<double> dfdx_all(double u, int32_t n_result, int32_t* ierr = nullptr) override {
        std::vector<double> result(n_result);
        fitpack_closed_curve_c_curve_all_derivatives(as<fitpack_closed_curve_c>(), u, ierr, result.data(), n_result);
        return result;
    }


    int32_t comm_size() const override {
        return fitpack_closed_curve_c_comm_size(as<fitpack_closed_curve_c>());
    }

    /**
     * @brief comm_pack
     */
    void comm_pack(std::vector<double>& buffer) const override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_closed_curve_c_comm_pack(as<fitpack_closed_curve_c>(), n, buffer.data());
    }


    /**
     * @brief comm_expand
     */
    void comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_closed_curve_c_comm_expand(as<fitpack_closed_curve_c>(), n, buffer.data());
    }


    double mse() const override {
        return fitpack_closed_curve_c_mse(as<fitpack_closed_curve_c>());
    }

    int32_t core_comm_size() const override {
        return fitpack_closed_curve_c_core_comm_size(as<fitpack_closed_curve_c>());
    }

    /**
     * @brief core_comm_pack
     */
    void core_comm_pack(std::vector<double>& buffer) const override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_closed_curve_c_core_comm_pack(as<fitpack_closed_curve_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    void core_comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_closed_curve_c_core_comm_expand(as<fitpack_closed_curve_c>(), n, buffer.data());
    }


    void destroy_base() override {
        fitpack_closed_curve_c_destroy_base(as<fitpack_closed_curve_c>());
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

protected:
    explicit fpClosedCurve(NoAlloc tag) : fpParametricCurve(tag) {}

};

#endif /* FPCLOSEDCURVE_HPP_INCLUDED */
