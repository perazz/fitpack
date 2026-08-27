/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fpSurface.hpp (class fpSurface)
!> @brief Standalone C++ wrapper for fitpack_surface (no fortran-arrays dependency)
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

#ifndef FPSURFACE_HPP_INCLUDED
#define FPSURFACE_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_surface_c.h"
#include "fpFitter.hpp"
#include "fpCurve.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <variant>
#include <optional>
#include <array>

static_assert(sizeof(fitpack_surface_c) == sizeof(fitpack_fitter_c),
    "C descriptor layout mismatch: fitpack_surface_c vs fitpack_fitter_c");

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_surface
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_surface, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 * Extends fpFitter (Fortran: extends(fitpack_fitter))
 */
class fpSurface : public fpFitter {
public:
    // ===========================================================================================
    // Constructors and Destructor
    // ===========================================================================================

    /**
     * @brief Default constructor - allocates new fitpack_surface
     */
    fpSurface() : fpFitter(NoAlloc{}) {
        fitpack_surface_c_allocate(as<fitpack_surface_c>(), nullptr);
    }

    /**
     * @brief Destructor - deallocates if owned
     */
    ~fpSurface() override {
        fitpack_surface_c_destroy(as<fitpack_surface_c>(), nullptr);
    }

    /**
     * @brief Copy constructor - deep copy
     */
    fpSurface(const fpSurface& other) : fpFitter(NoAlloc{}) {
        *as<fitpack_surface_c>() = fitpack_surface_c_null;
        fitpack_surface_c_copy(as<fitpack_surface_c>(), other.as<fitpack_surface_c>(), false, nullptr);
    }

    /**
     * @brief Copy assignment - deep copy
     */
    fpSurface& operator=(const fpSurface& other) {
        if (this != &other) {
            fitpack_surface_c_destroy(as<fitpack_surface_c>(), nullptr);
            fitpack_surface_c_copy(as<fitpack_surface_c>(), other.as<fitpack_surface_c>(), false, nullptr);
        }
        return *this;
    }

    /**
     * @brief Move constructor - transfer ownership
     */
    fpSurface(fpSurface&& other) noexcept : fpFitter(NoAlloc{}) {
        *as<fitpack_surface_c>() = fitpack_surface_c_null;
        fitpack_surface_c_move_alloc(as<fitpack_surface_c>(), other.as<fitpack_surface_c>(), nullptr);
    }

    /**
     * @brief Move assignment - transfer ownership
     */
    fpSurface& operator=(fpSurface&& other) noexcept {
        if (this != &other) {
            fitpack_surface_c_destroy(as<fitpack_surface_c>(), nullptr);
            fitpack_surface_c_move_alloc(as<fitpack_surface_c>(), other.as<fitpack_surface_c>(), nullptr);
        }
        return *this;
    }

    /**
     * @brief Construct from existing C wrapper (takes ownership if move=true)
     */
    explicit fpSurface(fitpack_surface_c& c_wrapper, bool move = false) : fpFitter(NoAlloc{}) {
        if (move) {
            fitpack_surface_c_move_alloc(as<fitpack_surface_c>(), &c_wrapper, nullptr);
        } else {
            *as<fitpack_surface_c>() = fitpack_surface_c_null;
            fitpack_surface_c_copy(as<fitpack_surface_c>(), &c_wrapper, false, nullptr);
        }
    }

    /**
     * @brief Non-owning view ctor: bit-copy the C handle, mark
     * is_pointer=true. Used by parent classes' polymorphic accessors
     * to wrap a base-class handle as the correct concrete subtype.
     */
    explicit fpSurface(fitpack_surface_c& c_wrapper, ViewTag) : fpFitter(NoAlloc{}) {
        *as<fitpack_surface_c>() = c_wrapper;
        as<fitpack_surface_c>()->is_pointer = true;
    }

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    bool is_allocated() const override {
        return as<fitpack_surface_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    bool is_pointer() const override {
        return as<fitpack_surface_c>()->is_pointer;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_surface_c").
     */
    const char* c_type_name() const override {
        return fitpack_surface_c_c_type_name(*as<fitpack_surface_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    const char* cpp_type_name() const override {
        return "fpSurface";
    }

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_surface_c& c_handle() {
        return *as<fitpack_surface_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_surface_c& c_handle() const {
        return *as<fitpack_surface_c>();
    }

    /**
     * @brief Construct an empty (non-owning) wrapper, for use as a view target.
     *
     * The returned object has a null handle (cptr==nullptr, is_pointer==false)
     * and does NOT allocate a Fortran object. Intended as a fill target for
     * member-view accessors emitted on parent types — the parent populates
     * the handle with c_loc() into its own storage and sets is_pointer=true.
     */
    static fpSurface make_view() {
        return fpSurface(NoAlloc{});
    }

    /**
     * @brief Upcast to parent type (reference, no copy)
     */
    fpFitter& as_parent() {
        return static_cast<fpFitter&>(*this);
    }
    const fpFitter& as_parent() const {
        return static_cast<const fpFitter&>(*this);
    }

    // ===========================================================================================
    // Method Wrappers (standalone — no fxArray dependency)
    // ===========================================================================================

    /**
     * @brief new_points
     */
    virtual void new_points(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double* w = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        fitpack_surface_c_new_points(as<fitpack_surface_c>(), n, x.data(), y.data(), z.data(), w);
    }


    /**
     * @brief new_fit
     */
    virtual int32_t new_fit(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double* w = nullptr, double* smoothing = nullptr, int32_t* order = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        return fitpack_surface_c_new_fit(as<fitpack_surface_c>(), n, x.data(), y.data(), z.data(), w, smoothing, order);
    }


    virtual int32_t fit(double* smoothing = nullptr, int32_t* order = nullptr, bool* keep_knots = nullptr) {
        return fitpack_surface_c_fit(as<fitpack_surface_c>(), smoothing, order, keep_knots);
    }

    virtual int32_t interpolate(bool* reset_knots = nullptr) {
        return fitpack_surface_c_interpolate(as<fitpack_surface_c>(), reset_knots);
    }

    virtual int32_t least_squares(double* smoothing = nullptr, bool* reset_knots = nullptr) {
        return fitpack_surface_c_least_squares(as<fitpack_surface_c>(), smoothing, reset_knots);
    }

    virtual double eval(double x, double y, int32_t* ierr = nullptr) {
        return fitpack_surface_c_surface_eval_one(as<fitpack_surface_c>(), x, y, ierr);
    }

    /**
     * @brief eval
     */
    virtual std::vector<double> eval(std::vector<double>& x, std::vector<double>& y, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n);
        int32_t n_result = 0;
        fitpack_surface_c_surface_eval_many(as<fitpack_surface_c>(), n, x.data(), y.data(), ierr, result.data(), &n_result, n);
        result.resize(n_result);
        return result;
    }


    /**
     * @brief eval_ongrid
     */
    virtual std::vector<double> eval_ongrid(std::vector<double>& x, std::vector<double>& y, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n * n);
        int32_t n_result = 0;
        fitpack_surface_c_surface_eval_gridded(as<fitpack_surface_c>(), n, x.data(), y.data(), ierr, result.data(), &n_result, n * n);
        result.resize(n_result);
        return result;
    }


    virtual double dfdx(double x, double y, int32_t dx, int32_t dy, int32_t* ierr = nullptr) {
        return fitpack_surface_c_surface_derivatives_one(as<fitpack_surface_c>(), x, y, dx, dy, ierr);
    }

    /**
     * @brief dfdx
     */
    virtual std::vector<double> dfdx(std::vector<double>& x, std::vector<double>& y, int32_t dx, int32_t dy, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n);
        int32_t n_result = 0;
        fitpack_surface_c_surface_derivatives_many(as<fitpack_surface_c>(), n, x.data(), y.data(), dx, dy, ierr, result.data(), &n_result, n);
        result.resize(n_result);
        return result;
    }


    /**
     * @brief dfdx_ongrid
     */
    virtual std::vector<double> dfdx_ongrid(std::vector<double>& x, std::vector<double>& y, int32_t dx, int32_t dy, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n * n);
        int32_t n_result = 0;
        fitpack_surface_c_surface_derivatives_gridded(as<fitpack_surface_c>(), n, x.data(), y.data(), dx, dy, ierr, result.data(), &n_result, n * n);
        result.resize(n_result);
        return result;
    }


    /**
     * @brief integral
     */
    virtual double integral(std::vector<double>& lower, std::vector<double>& upper) const {
        int32_t n = static_cast<int32_t>(lower.size());
        return fitpack_surface_c_integral(as<fitpack_surface_c>(), n, lower.data(), upper.data());
    }


    virtual fpCurve cross_section(double u, bool along_y, int32_t* ierr = nullptr) {
        fitpack_curve_c result_c;
        fitpack_surface_c_cross_section(as<fitpack_surface_c>(), u, along_y, ierr, &result_c);
        return fpCurve(result_c, true);
    }

    virtual fpSurface derivative_spline(int32_t nux, int32_t nuy, int32_t* ierr = nullptr) {
        fitpack_surface_c result_c;
        fitpack_surface_c_derivative_spline(as<fitpack_surface_c>(), nux, nuy, ierr, &result_c);
        return fpSurface(result_c, true);
    }

    int32_t comm_size() const override {
        return fitpack_surface_c_comm_size(as<fitpack_surface_c>());
    }

    /**
     * @brief comm_pack
     */
    void comm_pack(std::vector<double>& buffer) const override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_surface_c_comm_pack(as<fitpack_surface_c>(), n, buffer.data());
    }


    /**
     * @brief comm_expand
     */
    void comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_surface_c_comm_expand(as<fitpack_surface_c>(), n, buffer.data());
    }


    double mse() const override {
        return fitpack_surface_c_mse(as<fitpack_surface_c>());
    }

    int32_t core_comm_size() const override {
        return fitpack_surface_c_core_comm_size(as<fitpack_surface_c>());
    }

    /**
     * @brief core_comm_pack
     */
    void core_comm_pack(std::vector<double>& buffer) const override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_surface_c_core_comm_pack(as<fitpack_surface_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    void core_comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_surface_c_core_comm_expand(as<fitpack_surface_c>(), n, buffer.data());
    }


    void destroy_base() override {
        fitpack_surface_c_destroy_base(as<fitpack_surface_c>());
    }

    // ===========================================================================================
    // Component Array Accessors
    // ===========================================================================================

    /**
     * @brief Deep copy of component 'x' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> x_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_surface_c_getcomp_x(as<fitpack_surface_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'x'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through x()(i) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> x() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_surface_c_getcomp_x(as<fitpack_surface_c>(), &raw, extents);
        FX_SIZE bounds[2];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 1; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "x", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(1), bounds);
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    /**
     * @brief Deep copy of component 'y' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> y_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_surface_c_getcomp_y(as<fitpack_surface_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'y'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through y()(i) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> y() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_surface_c_getcomp_y(as<fitpack_surface_c>(), &raw, extents);
        FX_SIZE bounds[2];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 1; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "y", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(1), bounds);
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    /**
     * @brief Deep copy of component 'z' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> z_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_surface_c_getcomp_z(as<fitpack_surface_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'z'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through z()(i) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> z() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_surface_c_getcomp_z(as<fitpack_surface_c>(), &raw, extents);
        FX_SIZE bounds[2];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 1; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "z", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(1), bounds);
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    /**
     * @brief Deep copy of component 'w' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> w_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_surface_c_getcomp_w(as<fitpack_surface_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'w'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through w()(i) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> w() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_surface_c_getcomp_w(as<fitpack_surface_c>(), &raw, extents);
        FX_SIZE bounds[2];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 1; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "w", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(1), bounds);
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    /**
     * @brief Deep copy of component 'wrk2' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> wrk2_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_surface_c_getcomp_wrk2(as<fitpack_surface_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'wrk2'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through wrk2()(i) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> wrk2() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_surface_c_getcomp_wrk2(as<fitpack_surface_c>(), &raw, extents);
        FX_SIZE bounds[2];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 1; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "wrk2", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(1), bounds);
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    /**
     * @brief Deep copy of component 't' as a std::vector.
     *
     * The rank-2 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with t_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> t_vector() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_surface_c_getcomp_t(as<fitpack_surface_c>(), &raw, extents);
        const int64_t total = extents[0] * extents[1];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

    /**
     * @brief Extents of component 't', leading dimension first.
     * @return All zeros when the component is unallocated.
     */
    std::array<int64_t, 2> t_shape() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_surface_c_getcomp_t(as<fitpack_surface_c>(), &raw, extents);
        return { extents[0], extents[1] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 't'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through t()(i, j) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> t() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_surface_c_getcomp_t(as<fitpack_surface_c>(), &raw, extents);
        FX_SIZE bounds[4];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 2; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "t", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(2), bounds);
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    // ===========================================================================================
    // Scalar Property Accessors
    // ===========================================================================================

    int32_t& m() {
        return *fitpack_surface_c_ref_m(as<fitpack_surface_c>());
    }
    const int32_t& m() const {
        return *fitpack_surface_c_ref_m(as<fitpack_surface_c>());
    }

    int32_t& nmax() {
        return *fitpack_surface_c_ref_nmax(as<fitpack_surface_c>());
    }
    const int32_t& nmax() const {
        return *fitpack_surface_c_ref_nmax(as<fitpack_surface_c>());
    }

    int32_t& lwrk2() {
        return *fitpack_surface_c_ref_lwrk2(as<fitpack_surface_c>());
    }
    const int32_t& lwrk2() const {
        return *fitpack_surface_c_ref_lwrk2(as<fitpack_surface_c>());
    }

    int32_t& bc() {
        return *fitpack_surface_c_ref_bc(as<fitpack_surface_c>());
    }
    const int32_t& bc() const {
        return *fitpack_surface_c_ref_bc(as<fitpack_surface_c>());
    }

protected:
    explicit fpSurface(NoAlloc tag) : fpFitter(tag) {}

};

#endif /* FPSURFACE_HPP_INCLUDED */
