/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fpPeriodicCurve.hpp                                                                     **
 *   ** @brief Standalone C++ wrapper for fitpack_periodic_curve (no fortran-arrays dependency)
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FPPERIODICCURVE_HPP_INCLUDED
#define FPPERIODICCURVE_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_periodic_curve_c.h"
#include "fpCurve.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <variant>
#include <optional>
#include <cstring>

static_assert(sizeof(fitpack_periodic_curve_c) == sizeof(fitpack_curve_c),
    "C descriptor layout mismatch: fitpack_periodic_curve_c vs fitpack_curve_c");

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_periodic_curve
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_periodic_curve, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 * Extends fpCurve (Fortran: extends(fitpack_curve))
 */
class fpPeriodicCurve : public fpCurve {
public:
    // ===========================================================================================
    // Constructors and Destructor
    // ===========================================================================================

    /**
     * @brief Default constructor - allocates new fitpack_periodic_curve
     */
    fpPeriodicCurve() : fpCurve(NoAlloc{}) {
        fitpack_periodic_curve_c_allocate(as<fitpack_periodic_curve_c>(), nullptr);
    }

    /**
     * @brief Destructor - deallocates if owned
     */
    ~fpPeriodicCurve() override {
        fitpack_periodic_curve_c_destroy(as<fitpack_periodic_curve_c>(), nullptr);
    }

    /**
     * @brief Copy constructor - deep copy
     */
    fpPeriodicCurve(const fpPeriodicCurve& other) : fpCurve(NoAlloc{}) {
        *as<fitpack_periodic_curve_c>() = fitpack_periodic_curve_c_null;
        fitpack_periodic_curve_c_copy(as<fitpack_periodic_curve_c>(), other.as<fitpack_periodic_curve_c>(), false, nullptr);
    }

    /**
     * @brief Copy assignment - deep copy
     */
    fpPeriodicCurve& operator=(const fpPeriodicCurve& other) {
        if (this != &other) {
            fitpack_periodic_curve_c_destroy(as<fitpack_periodic_curve_c>(), nullptr);
            fitpack_periodic_curve_c_copy(as<fitpack_periodic_curve_c>(), other.as<fitpack_periodic_curve_c>(), false, nullptr);
        }
        return *this;
    }

    /**
     * @brief Move constructor - transfer ownership
     */
    fpPeriodicCurve(fpPeriodicCurve&& other) noexcept : fpCurve(NoAlloc{}) {
        *as<fitpack_periodic_curve_c>() = fitpack_periodic_curve_c_null;
        fitpack_periodic_curve_c_move_alloc(as<fitpack_periodic_curve_c>(), other.as<fitpack_periodic_curve_c>(), nullptr);
    }

    /**
     * @brief Move assignment - transfer ownership
     */
    fpPeriodicCurve& operator=(fpPeriodicCurve&& other) noexcept {
        if (this != &other) {
            fitpack_periodic_curve_c_destroy(as<fitpack_periodic_curve_c>(), nullptr);
            fitpack_periodic_curve_c_move_alloc(as<fitpack_periodic_curve_c>(), other.as<fitpack_periodic_curve_c>(), nullptr);
        }
        return *this;
    }

    /**
     * @brief Construct from existing C wrapper (takes ownership if move=true)
     */
    explicit fpPeriodicCurve(fitpack_periodic_curve_c& c_wrapper, bool move = false) : fpCurve(NoAlloc{}) {
        if (move) {
            fitpack_periodic_curve_c_move_alloc(as<fitpack_periodic_curve_c>(), &c_wrapper, nullptr);
        } else {
            *as<fitpack_periodic_curve_c>() = fitpack_periodic_curve_c_null;
            fitpack_periodic_curve_c_copy(as<fitpack_periodic_curve_c>(), &c_wrapper, false, nullptr);
        }
    }

    /**
     * @brief Non-owning view ctor: bit-copy the C handle, mark
     * is_pointer=true. Used by parent classes' polymorphic accessors
     * to wrap a base-class handle as the correct concrete subtype.
     */
    explicit fpPeriodicCurve(fitpack_periodic_curve_c& c_wrapper, ViewTag) : fpCurve(NoAlloc{}) {
        *as<fitpack_periodic_curve_c>() = c_wrapper;
        as<fitpack_periodic_curve_c>()->is_pointer = true;
    }

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    bool is_allocated() const override {
        return as<fitpack_periodic_curve_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    bool is_pointer() const override {
        return as<fitpack_periodic_curve_c>()->is_pointer;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_periodic_curve_c").
     */
    const char* c_type_name() const override {
        return fitpack_periodic_curve_c_c_type_name(*as<fitpack_periodic_curve_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    const char* cpp_type_name() const override {
        return "fpPeriodicCurve";
    }

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_periodic_curve_c& c_handle() {
        return *as<fitpack_periodic_curve_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_periodic_curve_c& c_handle() const {
        return *as<fitpack_periodic_curve_c>();
    }

    /**
     * @brief Construct an empty (non-owning) wrapper, for use as a view target.
     *
     * The returned object has a null handle (cptr==nullptr, is_pointer==false)
     * and does NOT allocate a Fortran object. Intended as a fill target for
     * member-view accessors emitted on parent types — the parent populates
     * the handle with c_loc() into its own storage and sets is_pointer=true.
     */
    static fpPeriodicCurve make_view() {
        return fpPeriodicCurve(NoAlloc{});
    }

    /**
     * @brief Upcast to parent type (reference, no copy)
     */
    fpCurve& as_parent() {
        return static_cast<fpCurve&>(*this);
    }
    const fpCurve& as_parent() const {
        return static_cast<const fpCurve&>(*this);
    }

    // ===========================================================================================
    // Method Wrappers (standalone — no fxArray dependency)
    // ===========================================================================================

    /**
     * @brief new_points
     */
    void new_points(std::vector<double>& x, std::vector<double>& y, double* w = nullptr) override {
        int32_t n = static_cast<int32_t>(x.size());
        fitpack_periodic_curve_c_new_points(as<fitpack_periodic_curve_c>(), n, x.data(), y.data(), w);
    }


    /**
     * @brief new_fit
     */
    int32_t new_fit(std::vector<double>& x, std::vector<double>& y, double* w = nullptr, double* smoothing = nullptr, int32_t* order = nullptr) const override {
        int32_t n = static_cast<int32_t>(x.size());
        return fitpack_periodic_curve_c_new_fit(as<fitpack_periodic_curve_c>(), n, x.data(), y.data(), w, smoothing, order);
    }


    int32_t fit(double* smoothing = nullptr, int32_t* order = nullptr, bool* keep_knots = nullptr) const override {
        return fitpack_periodic_curve_c_fit(as<fitpack_periodic_curve_c>(), smoothing, order, keep_knots);
    }

    int32_t interpolate(int32_t* order = nullptr, bool* reset_knots = nullptr) const override {
        return fitpack_periodic_curve_c_interpolate(as<fitpack_periodic_curve_c>(), order, reset_knots);
    }

    int32_t least_squares(double* smoothing = nullptr, bool* reset_knots = nullptr) const override {
        return fitpack_periodic_curve_c_least_squares(as<fitpack_periodic_curve_c>(), smoothing, reset_knots);
    }

    double eval(double x, int32_t& ierr) const override {
        return fitpack_periodic_curve_c_curve_eval_one(as<fitpack_periodic_curve_c>(), x, &ierr);
    }

    double eval(double x) const override {
        return fitpack_periodic_curve_c_curve_eval_one_noerr(as<fitpack_periodic_curve_c>(), x);
    }

    /**
     * @brief eval
     */
    std::vector<double> eval(std::vector<double>& x, int32_t& ierr) override {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n);
        int32_t n_result = 0;
        fitpack_periodic_curve_c_curve_eval_many(as<fitpack_periodic_curve_c>(), n, x.data(), &ierr, result.data(), &n_result, n);
        result.resize(n_result);
        return result;
    }


    /**
     * @brief eval
     */
    std::vector<double> eval(std::vector<double>& x) override {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n);
        int32_t n_result = 0;
        fitpack_periodic_curve_c_curve_eval_many_pure(as<fitpack_periodic_curve_c>(), n, x.data(), result.data(), &n_result, n);
        result.resize(n_result);
        return result;
    }


    double integral(double from, double to) const override {
        return fitpack_periodic_curve_c_integral(as<fitpack_periodic_curve_c>(), from, to);
    }

    /**
     * @brief fourier_coefficients
     */
    void fourier_coefficients(std::vector<double>& alpha, std::vector<double>& a, std::vector<double>& b, int32_t* ierr = nullptr) override {
        int32_t n = static_cast<int32_t>(alpha.size());
        fitpack_periodic_curve_c_fourier_coefficients(as<fitpack_periodic_curve_c>(), n, alpha.data(), a.data(), b.data(), ierr);
    }


    /**
     * @brief zeros
     */
    std::vector<double> zeros(int32_t max_size, int32_t* ierr = nullptr) override {
        std::vector<double> result(max_size);
        int32_t n_result = 0;
        fitpack_periodic_curve_c_zeros(as<fitpack_periodic_curve_c>(), ierr, result.data(), &n_result, max_size);
        result.resize(n_result);
        return result;
    }


    double dfdx(double x, int32_t order, int32_t* ierr = nullptr) const override {
        return fitpack_periodic_curve_c_curve_derivative(as<fitpack_periodic_curve_c>(), x, order, ierr);
    }

    /**
     * @brief dfdx
     */
    std::vector<double> dfdx(std::vector<double>& x, int32_t order, int32_t* ierr = nullptr) override {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n);
        int32_t n_result = 0;
        fitpack_periodic_curve_c_curve_derivatives(as<fitpack_periodic_curve_c>(), n, x.data(), order, ierr, result.data(), &n_result, n);
        result.resize(n_result);
        return result;
    }


    /**
     * @brief dfdx_all
     */
    std::vector<double> dfdx_all(double x, int32_t& ierr, int32_t n_result) override {
        std::vector<double> result(n_result);
        fitpack_periodic_curve_c_curve_all_derivatives(as<fitpack_periodic_curve_c>(), x, &ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx_all
     */
    std::vector<double> dfdx_all(double x, int32_t n_result) override {
        std::vector<double> result(n_result);
        fitpack_periodic_curve_c_curve_all_derivatives_pure(as<fitpack_periodic_curve_c>(), x, result.data(), n_result);
        return result;
    }


    void insert_knot(double x, int32_t* ierr = nullptr) override {
        fitpack_periodic_curve_c_curve_insert_knot_one(as<fitpack_periodic_curve_c>(), x, ierr);
    }

    /**
     * @brief insert_knot
     */
    void insert_knot(std::vector<double>& x, int32_t* ierr = nullptr) override {
        int32_t n = static_cast<int32_t>(x.size());
        fitpack_periodic_curve_c_curve_insert_knot_many(as<fitpack_periodic_curve_c>(), n, x.data(), ierr);
    }


    int32_t comm_size() const override {
        return fitpack_periodic_curve_c_comm_size(as<fitpack_periodic_curve_c>());
    }

    /**
     * @brief comm_pack
     */
    void comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_periodic_curve_c_comm_pack(as<fitpack_periodic_curve_c>(), n, buffer.data());
    }


    /**
     * @brief comm_expand
     */
    void comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_periodic_curve_c_comm_expand(as<fitpack_periodic_curve_c>(), n, buffer.data());
    }


    double mse() const override {
        return fitpack_periodic_curve_c_mse(as<fitpack_periodic_curve_c>());
    }

    int32_t core_comm_size() const override {
        return fitpack_periodic_curve_c_core_comm_size(as<fitpack_periodic_curve_c>());
    }

    /**
     * @brief core_comm_pack
     */
    void core_comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_periodic_curve_c_core_comm_pack(as<fitpack_periodic_curve_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    void core_comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_periodic_curve_c_core_comm_expand(as<fitpack_periodic_curve_c>(), n, buffer.data());
    }


    void destroy_base() override {
        fitpack_periodic_curve_c_destroy_base(as<fitpack_periodic_curve_c>());
    }

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

protected:
    explicit fpPeriodicCurve(NoAlloc tag) : fpCurve(tag) {}

};

#endif /* FPPERIODICCURVE_HPP_INCLUDED */
