/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fpGridSurface.hpp (class fpGridSurface)
!> @brief Standalone C++ wrapper for fitpack_grid_surface (no fortran-arrays dependency)
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

#ifndef FPGRIDSURFACE_HPP_INCLUDED
#define FPGRIDSURFACE_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_grid_surface_c.h"
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

static_assert(sizeof(fitpack_grid_surface_c) == sizeof(fitpack_fitter_c),
    "C descriptor layout mismatch: fitpack_grid_surface_c vs fitpack_fitter_c");

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_grid_surface
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_grid_surface, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 * Extends fpFitter (Fortran: extends(fitpack_fitter))
 */
class fpGridSurface : public fpFitter {
public:
    // ===========================================================================================
    // Constructors and Destructor
    // ===========================================================================================

    /**
     * @brief Default constructor - allocates new fitpack_grid_surface
     */
    fpGridSurface() : fpFitter(NoAlloc{}) {
        fitpack_grid_surface_c_allocate(as<fitpack_grid_surface_c>(), nullptr);
    }

    /**
     * @brief Destructor - deallocates if owned
     */
    ~fpGridSurface() override {
        fitpack_grid_surface_c_destroy(as<fitpack_grid_surface_c>(), nullptr);
    }

    /**
     * @brief Copy constructor - deep copy
     */
    fpGridSurface(const fpGridSurface& other) : fpFitter(NoAlloc{}) {
        *as<fitpack_grid_surface_c>() = fitpack_grid_surface_c_null;
        fitpack_grid_surface_c_copy(as<fitpack_grid_surface_c>(), other.as<fitpack_grid_surface_c>(), false, nullptr);
    }

    /**
     * @brief Copy assignment - deep copy
     */
    fpGridSurface& operator=(const fpGridSurface& other) {
        if (this != &other) {
            fitpack_grid_surface_c_destroy(as<fitpack_grid_surface_c>(), nullptr);
            fitpack_grid_surface_c_copy(as<fitpack_grid_surface_c>(), other.as<fitpack_grid_surface_c>(), false, nullptr);
        }
        return *this;
    }

    /**
     * @brief Move constructor - transfer ownership
     */
    fpGridSurface(fpGridSurface&& other) noexcept : fpFitter(NoAlloc{}) {
        *as<fitpack_grid_surface_c>() = fitpack_grid_surface_c_null;
        fitpack_grid_surface_c_move_alloc(as<fitpack_grid_surface_c>(), other.as<fitpack_grid_surface_c>(), nullptr);
    }

    /**
     * @brief Move assignment - transfer ownership
     */
    fpGridSurface& operator=(fpGridSurface&& other) noexcept {
        if (this != &other) {
            fitpack_grid_surface_c_destroy(as<fitpack_grid_surface_c>(), nullptr);
            fitpack_grid_surface_c_move_alloc(as<fitpack_grid_surface_c>(), other.as<fitpack_grid_surface_c>(), nullptr);
        }
        return *this;
    }

    /**
     * @brief Construct from existing C wrapper (takes ownership if move=true)
     */
    explicit fpGridSurface(fitpack_grid_surface_c& c_wrapper, bool move = false) : fpFitter(NoAlloc{}) {
        if (move) {
            fitpack_grid_surface_c_move_alloc(as<fitpack_grid_surface_c>(), &c_wrapper, nullptr);
        } else {
            *as<fitpack_grid_surface_c>() = fitpack_grid_surface_c_null;
            fitpack_grid_surface_c_copy(as<fitpack_grid_surface_c>(), &c_wrapper, false, nullptr);
        }
    }

    /**
     * @brief Non-owning view ctor: bit-copy the C handle, mark
     * is_pointer=true. Used by parent classes' polymorphic accessors
     * to wrap a base-class handle as the correct concrete subtype.
     */
    explicit fpGridSurface(fitpack_grid_surface_c& c_wrapper, ViewTag) : fpFitter(NoAlloc{}) {
        *as<fitpack_grid_surface_c>() = c_wrapper;
        as<fitpack_grid_surface_c>()->is_pointer = true;
    }

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    bool is_allocated() const override {
        return as<fitpack_grid_surface_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    bool is_pointer() const override {
        return as<fitpack_grid_surface_c>()->is_pointer;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_grid_surface_c").
     */
    const char* c_type_name() const override {
        return fitpack_grid_surface_c_c_type_name(*as<fitpack_grid_surface_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    const char* cpp_type_name() const override {
        return "fpGridSurface";
    }

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_grid_surface_c& c_handle() {
        return *as<fitpack_grid_surface_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_grid_surface_c& c_handle() const {
        return *as<fitpack_grid_surface_c>();
    }

    /**
     * @brief Construct an empty (non-owning) wrapper, for use as a view target.
     *
     * The returned object has a null handle (cptr==nullptr, is_pointer==false)
     * and does NOT allocate a Fortran object. Intended as a fill target for
     * member-view accessors emitted on parent types — the parent populates
     * the handle with c_loc() into its own storage and sets is_pointer=true.
     */
    static fpGridSurface make_view() {
        return fpGridSurface(NoAlloc{});
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
    virtual void new_points(std::vector<double>& x, std::vector<double>& y, int32_t z_n1, int32_t z_n2, std::vector<double>& z) {
        int32_t n = static_cast<int32_t>(x.size());
        fitpack_grid_surface_c_new_points(as<fitpack_grid_surface_c>(), n, x.data(), y.data(), z_n1, z_n2, z.data());
    }


    /**
     * @brief new_fit
     */
    virtual int32_t new_fit(std::vector<double>& x, std::vector<double>& y, int32_t z_n1, int32_t z_n2, std::vector<double>& z, double* smoothing = nullptr, int32_t* order = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        return fitpack_grid_surface_c_new_fit(as<fitpack_grid_surface_c>(), n, x.data(), y.data(), z_n1, z_n2, z.data(), smoothing, order);
    }


    virtual int32_t fit(double* smoothing = nullptr, int32_t* order = nullptr, bool* keep_knots = nullptr) {
        return fitpack_grid_surface_c_fit(as<fitpack_grid_surface_c>(), smoothing, order, keep_knots);
    }

    virtual int32_t least_squares(double* smoothing = nullptr, bool* reset_knots = nullptr) {
        return fitpack_grid_surface_c_least_squares(as<fitpack_grid_surface_c>(), smoothing, reset_knots);
    }

    virtual int32_t interpolate(bool* reset_knots = nullptr) {
        return fitpack_grid_surface_c_interpolate(as<fitpack_grid_surface_c>(), reset_knots);
    }

    virtual double eval(double x, double y, int32_t* ierr = nullptr) {
        return fitpack_grid_surface_c_gridsurf_eval_one(as<fitpack_grid_surface_c>(), x, y, ierr);
    }

    /**
     * @brief eval
     */
    virtual std::vector<double> eval(std::vector<double>& x, std::vector<double>& y, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n);
        int32_t n_result = 0;
        fitpack_grid_surface_c_gridsurf_eval_many(as<fitpack_grid_surface_c>(), n, x.data(), y.data(), ierr, result.data(), &n_result, n);
        result.resize(n_result);
        return result;
    }


    virtual double eval_ongrid(double x, double y, int32_t* ierr = nullptr) {
        return fitpack_grid_surface_c_gridded_eval_one(as<fitpack_grid_surface_c>(), x, y, ierr);
    }

    /**
     * @brief eval_ongrid
     */
    virtual std::vector<double> eval_ongrid(std::vector<double>& x, std::vector<double>& y, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n * n);
        int32_t n_result = 0;
        fitpack_grid_surface_c_gridded_eval_many(as<fitpack_grid_surface_c>(), n, x.data(), y.data(), ierr, result.data(), &n_result, n * n);
        result.resize(n_result);
        return result;
    }


    virtual double dfdx(double x, double y, int32_t dx, int32_t dy, int32_t* ierr = nullptr) {
        return fitpack_grid_surface_c_gridded_derivatives_one(as<fitpack_grid_surface_c>(), x, y, dx, dy, ierr);
    }

    /**
     * @brief dfdx
     */
    virtual std::vector<double> dfdx(std::vector<double>& x, std::vector<double>& y, int32_t dx, int32_t dy, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n);
        int32_t n_result = 0;
        fitpack_grid_surface_c_gridded_derivatives_many(as<fitpack_grid_surface_c>(), n, x.data(), y.data(), dx, dy, ierr, result.data(), &n_result, n);
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
        fitpack_grid_surface_c_gridded_derivatives_gridded(as<fitpack_grid_surface_c>(), n, x.data(), y.data(), dx, dy, ierr, result.data(), &n_result, n * n);
        result.resize(n_result);
        return result;
    }


    /**
     * @brief integral
     */
    virtual double integral(std::vector<double>& lower, std::vector<double>& upper) const {
        int32_t n = static_cast<int32_t>(lower.size());
        return fitpack_grid_surface_c_integral(as<fitpack_grid_surface_c>(), n, lower.data(), upper.data());
    }


    virtual fpCurve cross_section(double u, bool along_y, int32_t* ierr = nullptr) {
        fitpack_curve_c result_c;
        fitpack_grid_surface_c_cross_section(as<fitpack_grid_surface_c>(), u, along_y, ierr, &result_c);
        return fpCurve(result_c, true);
    }

    virtual fpGridSurface derivative_spline(int32_t nux, int32_t nuy, int32_t* ierr = nullptr) {
        fitpack_grid_surface_c result_c;
        fitpack_grid_surface_c_derivative_spline(as<fitpack_grid_surface_c>(), nux, nuy, ierr, &result_c);
        return fpGridSurface(result_c, true);
    }

    int32_t comm_size() const override {
        return fitpack_grid_surface_c_comm_size(as<fitpack_grid_surface_c>());
    }

    /**
     * @brief comm_pack
     */
    void comm_pack(std::vector<double>& buffer) const override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_grid_surface_c_comm_pack(as<fitpack_grid_surface_c>(), n, buffer.data());
    }


    /**
     * @brief comm_expand
     */
    void comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_grid_surface_c_comm_expand(as<fitpack_grid_surface_c>(), n, buffer.data());
    }


    double mse() const override {
        return fitpack_grid_surface_c_mse(as<fitpack_grid_surface_c>());
    }

    int32_t core_comm_size() const override {
        return fitpack_grid_surface_c_core_comm_size(as<fitpack_grid_surface_c>());
    }

    /**
     * @brief core_comm_pack
     */
    void core_comm_pack(std::vector<double>& buffer) const override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_grid_surface_c_core_comm_pack(as<fitpack_grid_surface_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    void core_comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_grid_surface_c_core_comm_expand(as<fitpack_grid_surface_c>(), n, buffer.data());
    }


    void destroy_base() override {
        fitpack_grid_surface_c_destroy_base(as<fitpack_grid_surface_c>());
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
        fitpack_grid_surface_c_getcomp_x(as<fitpack_grid_surface_c>(), &raw, extents);
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
        fitpack_grid_surface_c_getcomp_x(as<fitpack_grid_surface_c>(), &raw, extents);
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
        fitpack_grid_surface_c_getcomp_y(as<fitpack_grid_surface_c>(), &raw, extents);
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
        fitpack_grid_surface_c_getcomp_y(as<fitpack_grid_surface_c>(), &raw, extents);
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
     *
     * The rank-2 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with z_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> z_vector() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_grid_surface_c_getcomp_z(as<fitpack_grid_surface_c>(), &raw, extents);
        const int64_t total = extents[0] * extents[1];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

    /**
     * @brief Extents of component 'z', leading dimension first.
     * @return All zeros when the component is unallocated.
     */
    std::array<int64_t, 2> z_shape() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_grid_surface_c_getcomp_z(as<fitpack_grid_surface_c>(), &raw, extents);
        return { extents[0], extents[1] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'z'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through z()(i, j) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> z() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_grid_surface_c_getcomp_z(as<fitpack_grid_surface_c>(), &raw, extents);
        FX_SIZE bounds[4];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 2; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "z", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(2), bounds);
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
        fitpack_grid_surface_c_getcomp_t(as<fitpack_grid_surface_c>(), &raw, extents);
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
        fitpack_grid_surface_c_getcomp_t(as<fitpack_grid_surface_c>(), &raw, extents);
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
        fitpack_grid_surface_c_getcomp_t(as<fitpack_grid_surface_c>(), &raw, extents);
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

    int32_t& nmax() {
        return *fitpack_grid_surface_c_ref_nmax(as<fitpack_grid_surface_c>());
    }
    const int32_t& nmax() const {
        return *fitpack_grid_surface_c_ref_nmax(as<fitpack_grid_surface_c>());
    }

protected:
    explicit fpGridSurface(NoAlloc tag) : fpFitter(tag) {}

};

#endif /* FPGRIDSURFACE_HPP_INCLUDED */
