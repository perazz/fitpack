/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fpGridSpline.hpp (class fpGridSpline)
!> @brief Standalone C++ wrapper for fitpack_gridded_spline (no fortran-arrays dependency)
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

#ifndef FPGRIDSPLINE_HPP_INCLUDED
#define FPGRIDSPLINE_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_gridded_spline_c.h"
#include "fpFitter.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <variant>
#include <optional>
#include <array>
#include "fpPoint.hpp"

static_assert(sizeof(fitpack_gridded_spline_c) == sizeof(fitpack_fitter_c),
    "C descriptor layout mismatch: fitpack_gridded_spline_c vs fitpack_fitter_c");

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_gridded_spline
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_gridded_spline, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 * Extends fpFitter (Fortran: extends(fitpack_fitter))
 */
class fpGridSpline : public fpFitter {
public:
    // ===========================================================================================
    // Constructors and Destructor
    // ===========================================================================================

    /**
     * @brief Default constructor - allocates new fitpack_gridded_spline
     */
    fpGridSpline() : fpFitter(NoAlloc{}) {
        fitpack_gridded_spline_c_allocate(as<fitpack_gridded_spline_c>(), nullptr);
    }

    /**
     * @brief Destructor - deallocates if owned
     */
    ~fpGridSpline() override {
        fitpack_gridded_spline_c_destroy(as<fitpack_gridded_spline_c>(), nullptr);
    }

    /**
     * @brief Copy constructor - deep copy
     */
    fpGridSpline(const fpGridSpline& other) : fpFitter(NoAlloc{}) {
        *as<fitpack_gridded_spline_c>() = fitpack_gridded_spline_c_null;
        fitpack_gridded_spline_c_copy(as<fitpack_gridded_spline_c>(), other.as<fitpack_gridded_spline_c>(), false, nullptr);
    }

    /**
     * @brief Copy assignment - deep copy
     */
    fpGridSpline& operator=(const fpGridSpline& other) {
        if (this != &other) {
            fitpack_gridded_spline_c_destroy(as<fitpack_gridded_spline_c>(), nullptr);
            fitpack_gridded_spline_c_copy(as<fitpack_gridded_spline_c>(), other.as<fitpack_gridded_spline_c>(), false, nullptr);
        }
        return *this;
    }

    /**
     * @brief Move constructor - transfer ownership
     */
    fpGridSpline(fpGridSpline&& other) noexcept : fpFitter(NoAlloc{}) {
        *as<fitpack_gridded_spline_c>() = fitpack_gridded_spline_c_null;
        fitpack_gridded_spline_c_move_alloc(as<fitpack_gridded_spline_c>(), other.as<fitpack_gridded_spline_c>(), nullptr);
    }

    /**
     * @brief Move assignment - transfer ownership
     */
    fpGridSpline& operator=(fpGridSpline&& other) noexcept {
        if (this != &other) {
            fitpack_gridded_spline_c_destroy(as<fitpack_gridded_spline_c>(), nullptr);
            fitpack_gridded_spline_c_move_alloc(as<fitpack_gridded_spline_c>(), other.as<fitpack_gridded_spline_c>(), nullptr);
        }
        return *this;
    }

    /**
     * @brief Construct from existing C wrapper (takes ownership if move=true)
     */
    explicit fpGridSpline(fitpack_gridded_spline_c& c_wrapper, bool move = false) : fpFitter(NoAlloc{}) {
        if (move) {
            fitpack_gridded_spline_c_move_alloc(as<fitpack_gridded_spline_c>(), &c_wrapper, nullptr);
        } else {
            *as<fitpack_gridded_spline_c>() = fitpack_gridded_spline_c_null;
            fitpack_gridded_spline_c_copy(as<fitpack_gridded_spline_c>(), &c_wrapper, false, nullptr);
        }
    }

    /**
     * @brief Non-owning view ctor: bit-copy the C handle, mark
     * is_pointer=true. Used by parent classes' polymorphic accessors
     * to wrap a base-class handle as the correct concrete subtype.
     */
    explicit fpGridSpline(fitpack_gridded_spline_c& c_wrapper, ViewTag) : fpFitter(NoAlloc{}) {
        *as<fitpack_gridded_spline_c>() = c_wrapper;
        as<fitpack_gridded_spline_c>()->is_pointer = true;
    }

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    bool is_allocated() const override {
        return as<fitpack_gridded_spline_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    bool is_pointer() const override {
        return as<fitpack_gridded_spline_c>()->is_pointer;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_gridded_spline_c").
     */
    const char* c_type_name() const override {
        return fitpack_gridded_spline_c_c_type_name(*as<fitpack_gridded_spline_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    const char* cpp_type_name() const override {
        return "fpGridSpline";
    }

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_gridded_spline_c& c_handle() {
        return *as<fitpack_gridded_spline_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_gridded_spline_c& c_handle() const {
        return *as<fitpack_gridded_spline_c>();
    }

    /**
     * @brief Construct an empty (non-owning) wrapper, for use as a view target.
     *
     * The returned object has a null handle (cptr==nullptr, is_pointer==false)
     * and does NOT allocate a Fortran object. Intended as a fill target for
     * member-view accessors emitted on parent types — the parent populates
     * the handle with c_loc() into its own storage and sets is_pointer=true.
     */
    static fpGridSpline make_view() {
        return fpGridSpline(NoAlloc{});
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
    virtual void new_points(int32_t n, int32_t xg_n1, int32_t xg_n2, std::vector<double>& xg, std::vector<double>& z, int32_t* m = nullptr, bool* row_major = nullptr, int32_t* order = nullptr) {
        fitpack_gridded_spline_c_new_points(as<fitpack_gridded_spline_c>(), xg_n1, xg_n2, xg.data(), static_cast<int32_t>(z.size()), z.data(), n, m, row_major, order);
    }


    /**
     * @brief new_fit
     */
    virtual int32_t new_fit(int32_t n, int32_t xg_n1, int32_t xg_n2, std::vector<double>& xg, std::vector<double>& z, int32_t* m = nullptr, bool* row_major = nullptr, int32_t* order = nullptr, double* smoothing = nullptr) {
        return fitpack_gridded_spline_c_new_fit(as<fitpack_gridded_spline_c>(), xg_n1, xg_n2, xg.data(), static_cast<int32_t>(z.size()), z.data(), n, m, row_major, order, smoothing);
    }


    virtual int32_t fit(double* smoothing = nullptr, int32_t* order = nullptr) {
        return fitpack_gridded_spline_c_fit(as<fitpack_gridded_spline_c>(), smoothing, order);
    }

    virtual int32_t least_squares(double* smoothing = nullptr, bool* reset_knots = nullptr) {
        return fitpack_gridded_spline_c_least_squares(as<fitpack_gridded_spline_c>(), smoothing, reset_knots);
    }

    virtual int32_t interpolate() {
        return fitpack_gridded_spline_c_interpolate(as<fitpack_gridded_spline_c>());
    }

    /**
     * @brief eval_ongrid
     */
    virtual std::vector<double> eval_ongrid(int32_t xg_n1, int32_t xg_n2, std::vector<double>& xg, std::vector<int32_t>& m, int32_t n_result, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(m.size());
        std::vector<double> result(n_result);
        fitpack_gridded_spline_c_eval_ongrid(as<fitpack_gridded_spline_c>(), xg_n1, xg_n2, xg.data(), n, m.data(), ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief eval
     */
    virtual double eval(std::vector<double>& x, int32_t* ierr = nullptr) const {
        int32_t n = static_cast<int32_t>(x.size());
        return fitpack_gridded_spline_c_grid_eval_one(as<fitpack_gridded_spline_c>(), n, x.data(), ierr);
    }


    /**
     * @brief eval
     */
    virtual std::vector<double> eval(int32_t xp_n1, int32_t xp_n2, std::vector<double>& xp, int32_t* ierr = nullptr) {
        std::vector<double> result(xp_n2);
        int32_t n_result = 0;
        fitpack_gridded_spline_c_grid_eval_many(as<fitpack_gridded_spline_c>(), xp_n1, xp_n2, xp.data(), ierr, result.data(), &n_result, xp_n2);
        result.resize(n_result);
        return result;
    }


    /**
     * @brief dfdx_ongrid
     */
    virtual std::vector<double> dfdx_ongrid(int32_t xg_n1, int32_t xg_n2, std::vector<double>& xg, std::vector<int32_t>& m, std::vector<int32_t>& nu, int32_t n_result, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(m.size());
        std::vector<double> result(n_result);
        fitpack_gridded_spline_c_dfdx_ongrid(as<fitpack_gridded_spline_c>(), xg_n1, xg_n2, xg.data(), n, m.data(), nu.data(), ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx
     */
    virtual double dfdx(std::vector<double>& x, std::vector<int32_t>& nu, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        return fitpack_gridded_spline_c_grid_derivatives_one(as<fitpack_gridded_spline_c>(), n, x.data(), nu.data(), ierr);
    }


    /**
     * @brief dfdx
     */
    virtual std::vector<double> dfdx(int32_t xp_n1, int32_t xp_n2, std::vector<double>& xp, std::vector<int32_t>& nu, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(nu.size());
        std::vector<double> result(xp_n2);
        int32_t n_result = 0;
        fitpack_gridded_spline_c_grid_derivatives_many(as<fitpack_gridded_spline_c>(), xp_n1, xp_n2, xp.data(), n, nu.data(), ierr, result.data(), &n_result, xp_n2);
        result.resize(n_result);
        return result;
    }


    /**
     * @brief integral
     */
    virtual double integral(std::vector<double>& lower, std::vector<double>& upper) const {
        int32_t n = static_cast<int32_t>(lower.size());
        return fitpack_gridded_spline_c_integral(as<fitpack_gridded_spline_c>(), n, lower.data(), upper.data());
    }


    virtual fpGridSpline cross_section(int32_t ax, double u, int32_t* ierr = nullptr) {
        fitpack_gridded_spline_c result_c;
        fitpack_gridded_spline_c_cross_section(as<fitpack_gridded_spline_c>(), ax, u, ierr, &result_c);
        return fpGridSpline(result_c, true);
    }

    /**
     * @brief derivative_spline
     */
    virtual fpGridSpline derivative_spline(std::vector<int32_t>& nu, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(nu.size());
        fitpack_gridded_spline_c result_c;
        fitpack_gridded_spline_c_derivative_spline(as<fitpack_gridded_spline_c>(), n, nu.data(), ierr, &result_c);
        return fpGridSpline(result_c, true);
    }


    int32_t comm_size() const override {
        return fitpack_gridded_spline_c_comm_size(as<fitpack_gridded_spline_c>());
    }

    /**
     * @brief comm_pack
     */
    void comm_pack(std::vector<double>& buffer) const override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_gridded_spline_c_comm_pack(as<fitpack_gridded_spline_c>(), n, buffer.data());
    }


    /**
     * @brief comm_expand
     */
    void comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_gridded_spline_c_comm_expand(as<fitpack_gridded_spline_c>(), n, buffer.data());
    }


    double mse() const override {
        return fitpack_gridded_spline_c_mse(as<fitpack_gridded_spline_c>());
    }

    int32_t core_comm_size() const override {
        return fitpack_gridded_spline_c_core_comm_size(as<fitpack_gridded_spline_c>());
    }

    /**
     * @brief core_comm_pack
     */
    void core_comm_pack(std::vector<double>& buffer) const override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_gridded_spline_c_core_comm_pack(as<fitpack_gridded_spline_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    void core_comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_gridded_spline_c_core_comm_expand(as<fitpack_gridded_spline_c>(), n, buffer.data());
    }


    void destroy_base() override {
        fitpack_gridded_spline_c_destroy_base(as<fitpack_gridded_spline_c>());
    }

    // ===========================================================================================
    // Component Array Accessors
    // ===========================================================================================

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
        fitpack_gridded_spline_c_getcomp_t(as<fitpack_gridded_spline_c>(), &raw, extents);
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
        fitpack_gridded_spline_c_getcomp_t(as<fitpack_gridded_spline_c>(), &raw, extents);
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
        fitpack_gridded_spline_c_getcomp_t(as<fitpack_gridded_spline_c>(), &raw, extents);
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

    /**
     * @brief Deep copy of component 'xg' as a std::vector.
     *
     * The rank-2 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with xg_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> xg_vector() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_gridded_spline_c_getcomp_xg(as<fitpack_gridded_spline_c>(), &raw, extents);
        const int64_t total = extents[0] * extents[1];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

    /**
     * @brief Extents of component 'xg', leading dimension first.
     * @return All zeros when the component is unallocated.
     */
    std::array<int64_t, 2> xg_shape() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_gridded_spline_c_getcomp_xg(as<fitpack_gridded_spline_c>(), &raw, extents);
        return { extents[0], extents[1] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'xg'.
     *
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through xg()(i, j) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> xg() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_gridded_spline_c_getcomp_xg(as<fitpack_gridded_spline_c>(), &raw, extents);
        FX_SIZE bounds[4];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 2; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "xg", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(2), bounds);
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
        fitpack_gridded_spline_c_getcomp_z(as<fitpack_gridded_spline_c>(), &raw, extents);
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
        fitpack_gridded_spline_c_getcomp_z(as<fitpack_gridded_spline_c>(), &raw, extents);
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

    // ===========================================================================================
    // Scalar Property Accessors
    // ===========================================================================================

    int32_t& dims() {
        return *fitpack_gridded_spline_c_ref_dims(as<fitpack_gridded_spline_c>());
    }
    const int32_t& dims() const {
        return *fitpack_gridded_spline_c_ref_dims(as<fitpack_gridded_spline_c>());
    }

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
    //! @note Not const: the Fortran receiver of the raw dfdx is intent(inout).
    template <FP_SIZE dim>
    FP_REAL dfdx(const fpPoint<dim>& x, const std::vector<FP_SIZE>& nu, FP_FLAG* ierr = nullptr)
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

protected:
    explicit fpGridSpline(NoAlloc tag) : fpFitter(tag) {}

};

#endif /* FPGRIDSPLINE_HPP_INCLUDED */
