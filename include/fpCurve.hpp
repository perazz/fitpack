/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fpCurve.hpp                                                                     **
 *   ** @brief Standalone C++ wrapper for fitpack_curve (no fortran-arrays dependency)
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FPCURVE_HPP_INCLUDED
#define FPCURVE_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_curve_c.h"
#include "fxFitpackFitter.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <variant>
#include <optional>
#include <cstring>
#include <array>
#include "fpPoint.hpp"

static_assert(sizeof(fitpack_curve_c) == sizeof(fitpack_fitter_c),
    "C descriptor layout mismatch: fitpack_curve_c vs fitpack_fitter_c");

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_curve
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_curve, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 * Extends fxFitpackFitter (Fortran: extends(fitpack_fitter))
 */
class fpCurve : public fxFitpackFitter {
public:
    // ===========================================================================================
    // Constructors and Destructor
    // ===========================================================================================

    /**
     * @brief Default constructor - allocates new fitpack_curve
     */
    fpCurve() : fxFitpackFitter(NoAlloc{}) {
        fitpack_curve_c_allocate(as<fitpack_curve_c>(), nullptr);
    }

    /**
     * @brief Destructor - deallocates if owned
     */
    ~fpCurve() override {
        fitpack_curve_c_destroy(as<fitpack_curve_c>(), nullptr);
    }

    /**
     * @brief Copy constructor - deep copy
     */
    fpCurve(const fpCurve& other) : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_curve_c>() = fitpack_curve_c_null;
        fitpack_curve_c_copy(as<fitpack_curve_c>(), other.as<fitpack_curve_c>(), false, nullptr);
    }

    /**
     * @brief Copy assignment - deep copy
     */
    fpCurve& operator=(const fpCurve& other) {
        if (this != &other) {
            fitpack_curve_c_destroy(as<fitpack_curve_c>(), nullptr);
            fitpack_curve_c_copy(as<fitpack_curve_c>(), other.as<fitpack_curve_c>(), false, nullptr);
        }
        return *this;
    }

    /**
     * @brief Move constructor - transfer ownership
     */
    fpCurve(fpCurve&& other) noexcept : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_curve_c>() = fitpack_curve_c_null;
        fitpack_curve_c_move_alloc(as<fitpack_curve_c>(), other.as<fitpack_curve_c>(), nullptr);
    }

    /**
     * @brief Move assignment - transfer ownership
     */
    fpCurve& operator=(fpCurve&& other) noexcept {
        if (this != &other) {
            fitpack_curve_c_destroy(as<fitpack_curve_c>(), nullptr);
            fitpack_curve_c_move_alloc(as<fitpack_curve_c>(), other.as<fitpack_curve_c>(), nullptr);
        }
        return *this;
    }

    /**
     * @brief Construct from existing C wrapper (takes ownership if move=true)
     */
    explicit fpCurve(fitpack_curve_c& c_wrapper, bool move = false) : fxFitpackFitter(NoAlloc{}) {
        if (move) {
            fitpack_curve_c_move_alloc(as<fitpack_curve_c>(), &c_wrapper, nullptr);
        } else {
            *as<fitpack_curve_c>() = fitpack_curve_c_null;
            fitpack_curve_c_copy(as<fitpack_curve_c>(), &c_wrapper, false, nullptr);
        }
    }

    /**
     * @brief Non-owning view ctor: bit-copy the C handle, mark
     * is_pointer=true. Used by parent classes' polymorphic accessors
     * to wrap a base-class handle as the correct concrete subtype.
     */
    explicit fpCurve(fitpack_curve_c& c_wrapper, ViewTag) : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_curve_c>() = c_wrapper;
        as<fitpack_curve_c>()->is_pointer = true;
    }

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    bool is_allocated() const override {
        return as<fitpack_curve_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    bool is_pointer() const override {
        return as<fitpack_curve_c>()->is_pointer;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_curve_c").
     */
    const char* c_type_name() const override {
        return fitpack_curve_c_c_type_name(*as<fitpack_curve_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    const char* cpp_type_name() const override {
        return "fpCurve";
    }

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_curve_c& c_handle() {
        return *as<fitpack_curve_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_curve_c& c_handle() const {
        return *as<fitpack_curve_c>();
    }

    /**
     * @brief Construct an empty (non-owning) wrapper, for use as a view target.
     *
     * The returned object has a null handle (cptr==nullptr, is_pointer==false)
     * and does NOT allocate a Fortran object. Intended as a fill target for
     * member-view accessors emitted on parent types — the parent populates
     * the handle with c_loc() into its own storage and sets is_pointer=true.
     */
    static fpCurve make_view() {
        return fpCurve(NoAlloc{});
    }

    /**
     * @brief Upcast to parent type (reference, no copy)
     */
    fxFitpackFitter& as_parent() {
        return static_cast<fxFitpackFitter&>(*this);
    }
    const fxFitpackFitter& as_parent() const {
        return static_cast<const fxFitpackFitter&>(*this);
    }

    // ===========================================================================================
    // Method Wrappers (standalone — no fxArray dependency)
    // ===========================================================================================

    /**
     * @brief new_points
     */
    virtual void new_points(std::vector<double>& x, std::vector<double>& y, double* w = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        fitpack_curve_c_new_points(as<fitpack_curve_c>(), n, x.data(), y.data(), w);
    }


    /**
     * @brief new_fit
     */
    virtual int32_t new_fit(std::vector<double>& x, std::vector<double>& y, double* w = nullptr, double* smoothing = nullptr, int32_t* order = nullptr) const {
        int32_t n = static_cast<int32_t>(x.size());
        return fitpack_curve_c_new_fit(as<fitpack_curve_c>(), n, x.data(), y.data(), w, smoothing, order);
    }


    virtual int32_t fit(double* smoothing = nullptr, int32_t* order = nullptr, bool* keep_knots = nullptr) const {
        return fitpack_curve_c_fit(as<fitpack_curve_c>(), smoothing, order, keep_knots);
    }

    virtual int32_t interpolate(int32_t* order = nullptr, bool* reset_knots = nullptr) const {
        return fitpack_curve_c_interpolate(as<fitpack_curve_c>(), order, reset_knots);
    }

    virtual int32_t least_squares(double* smoothing = nullptr, bool* reset_knots = nullptr) const {
        return fitpack_curve_c_least_squares(as<fitpack_curve_c>(), smoothing, reset_knots);
    }

    virtual double eval(double x, int32_t& ierr) const {
        return fitpack_curve_c_curve_eval_one(as<fitpack_curve_c>(), x, &ierr);
    }

    virtual double eval(double x) const {
        return fitpack_curve_c_curve_eval_one_noerr(as<fitpack_curve_c>(), x);
    }

    /**
     * @brief eval
     */
    virtual std::vector<double> eval(std::vector<double>& x, int32_t& ierr) {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n);
        int32_t n_result = 0;
        fitpack_curve_c_curve_eval_many(as<fitpack_curve_c>(), n, x.data(), &ierr, result.data(), &n_result, n);
        result.resize(n_result);
        return result;
    }


    /**
     * @brief eval
     */
    virtual std::vector<double> eval(std::vector<double>& x) {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n);
        int32_t n_result = 0;
        fitpack_curve_c_curve_eval_many_pure(as<fitpack_curve_c>(), n, x.data(), result.data(), &n_result, n);
        result.resize(n_result);
        return result;
    }


    virtual double integral(double from, double to) const {
        return fitpack_curve_c_integral(as<fitpack_curve_c>(), from, to);
    }

    /**
     * @brief fourier_coefficients
     */
    virtual void fourier_coefficients(std::vector<double>& alpha, std::vector<double>& a, std::vector<double>& b, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(alpha.size());
        fitpack_curve_c_fourier_coefficients(as<fitpack_curve_c>(), n, alpha.data(), a.data(), b.data(), ierr);
    }


    /**
     * @brief zeros
     */
    virtual std::vector<double> zeros(int32_t max_size, int32_t* ierr = nullptr) {
        std::vector<double> result(max_size);
        int32_t n_result = 0;
        fitpack_curve_c_zeros(as<fitpack_curve_c>(), ierr, result.data(), &n_result, max_size);
        result.resize(n_result);
        return result;
    }


    virtual double dfdx(double x, int32_t order, int32_t* ierr = nullptr) const {
        return fitpack_curve_c_curve_derivative(as<fitpack_curve_c>(), x, order, ierr);
    }

    /**
     * @brief dfdx
     */
    virtual std::vector<double> dfdx(std::vector<double>& x, int32_t order, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        std::vector<double> result(n);
        int32_t n_result = 0;
        fitpack_curve_c_curve_derivatives(as<fitpack_curve_c>(), n, x.data(), order, ierr, result.data(), &n_result, n);
        result.resize(n_result);
        return result;
    }


    /**
     * @brief dfdx_all
     */
    virtual std::vector<double> dfdx_all(double x, int32_t& ierr, int32_t n_result) {
        std::vector<double> result(n_result);
        fitpack_curve_c_curve_all_derivatives(as<fitpack_curve_c>(), x, &ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx_all
     */
    virtual std::vector<double> dfdx_all(double x, int32_t n_result) {
        std::vector<double> result(n_result);
        fitpack_curve_c_curve_all_derivatives_pure(as<fitpack_curve_c>(), x, result.data(), n_result);
        return result;
    }


    virtual void insert_knot(double x, int32_t* ierr = nullptr) {
        fitpack_curve_c_curve_insert_knot_one(as<fitpack_curve_c>(), x, ierr);
    }

    /**
     * @brief insert_knot
     */
    virtual void insert_knot(std::vector<double>& x, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(x.size());
        fitpack_curve_c_curve_insert_knot_many(as<fitpack_curve_c>(), n, x.data(), ierr);
    }


    int32_t comm_size() const override {
        return fitpack_curve_c_comm_size(as<fitpack_curve_c>());
    }

    /**
     * @brief comm_pack
     */
    void comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_curve_c_comm_pack(as<fitpack_curve_c>(), n, buffer.data());
    }


    /**
     * @brief comm_expand
     */
    void comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_curve_c_comm_expand(as<fitpack_curve_c>(), n, buffer.data());
    }


    double mse() const override {
        return fitpack_curve_c_mse(as<fitpack_curve_c>());
    }

    int32_t core_comm_size() const override {
        return fitpack_curve_c_core_comm_size(as<fitpack_curve_c>());
    }

    /**
     * @brief core_comm_pack
     */
    void core_comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_curve_c_core_comm_pack(as<fitpack_curve_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    void core_comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_curve_c_core_comm_expand(as<fitpack_curve_c>(), n, buffer.data());
    }


    void destroy_base() override {
        fitpack_curve_c_destroy_base(as<fitpack_curve_c>());
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
        fitpack_curve_c_getcomp_x(as<fitpack_curve_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'x'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through x()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> x() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_curve_c_getcomp_x(as<fitpack_curve_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "x", FX_LEN_NAME - 1);
        descr.base_address = static_cast<void*>(raw);
        descr.type = getCFITypeFlag<double>();
        descr.elem_bytes = static_cast<FX_SIZE>(sizeof(double));
        descr.rank = static_cast<FX_RANK>(1);
        descr.is_pointer = true;   // non-owning: the object still owns the storage
        descr.is_slice = false;
        descr.attribute = static_cast<FX_ATTR>(FX_ATTR_POINTER);
        FX_SIZE stride = descr.elem_bytes;
        for (int k = 0; k < 1; ++k) {
            descr.dim[k].lower_bound = 0;   // C 0-based; Fortran lbound 1
            descr.dim[k].extent = static_cast<FX_SIZE>(extents[k]);
            descr.dim[k].stride_bytes = stride;
            stride *= descr.dim[k].extent;
        }
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
        fitpack_curve_c_getcomp_y(as<fitpack_curve_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'y'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through y()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> y() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_curve_c_getcomp_y(as<fitpack_curve_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "y", FX_LEN_NAME - 1);
        descr.base_address = static_cast<void*>(raw);
        descr.type = getCFITypeFlag<double>();
        descr.elem_bytes = static_cast<FX_SIZE>(sizeof(double));
        descr.rank = static_cast<FX_RANK>(1);
        descr.is_pointer = true;   // non-owning: the object still owns the storage
        descr.is_slice = false;
        descr.attribute = static_cast<FX_ATTR>(FX_ATTR_POINTER);
        FX_SIZE stride = descr.elem_bytes;
        for (int k = 0; k < 1; ++k) {
            descr.dim[k].lower_bound = 0;   // C 0-based; Fortran lbound 1
            descr.dim[k].extent = static_cast<FX_SIZE>(extents[k]);
            descr.dim[k].stride_bytes = stride;
            stride *= descr.dim[k].extent;
        }
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    /**
     * @brief Deep copy of component 'sp' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> sp_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_curve_c_getcomp_sp(as<fitpack_curve_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'sp'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through sp()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> sp() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_curve_c_getcomp_sp(as<fitpack_curve_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "sp", FX_LEN_NAME - 1);
        descr.base_address = static_cast<void*>(raw);
        descr.type = getCFITypeFlag<double>();
        descr.elem_bytes = static_cast<FX_SIZE>(sizeof(double));
        descr.rank = static_cast<FX_RANK>(1);
        descr.is_pointer = true;   // non-owning: the object still owns the storage
        descr.is_slice = false;
        descr.attribute = static_cast<FX_ATTR>(FX_ATTR_POINTER);
        FX_SIZE stride = descr.elem_bytes;
        for (int k = 0; k < 1; ++k) {
            descr.dim[k].lower_bound = 0;   // C 0-based; Fortran lbound 1
            descr.dim[k].extent = static_cast<FX_SIZE>(extents[k]);
            descr.dim[k].stride_bytes = stride;
            stride *= descr.dim[k].extent;
        }
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
        fitpack_curve_c_getcomp_w(as<fitpack_curve_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'w'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through w()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> w() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_curve_c_getcomp_w(as<fitpack_curve_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "w", FX_LEN_NAME - 1);
        descr.base_address = static_cast<void*>(raw);
        descr.type = getCFITypeFlag<double>();
        descr.elem_bytes = static_cast<FX_SIZE>(sizeof(double));
        descr.rank = static_cast<FX_RANK>(1);
        descr.is_pointer = true;   // non-owning: the object still owns the storage
        descr.is_slice = false;
        descr.attribute = static_cast<FX_ATTR>(FX_ATTR_POINTER);
        FX_SIZE stride = descr.elem_bytes;
        for (int k = 0; k < 1; ++k) {
            descr.dim[k].lower_bound = 0;   // C 0-based; Fortran lbound 1
            descr.dim[k].extent = static_cast<FX_SIZE>(extents[k]);
            descr.dim[k].stride_bytes = stride;
            stride *= descr.dim[k].extent;
        }
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    /**
     * @brief Deep copy of component 'wrk_fou' as a std::vector.
     *
     * The rank-2 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with wrk_fou_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> wrk_fou_vector() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_curve_c_getcomp_wrk_fou(as<fitpack_curve_c>(), &raw, extents);
        const int64_t total = extents[0] * extents[1];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

    /**
     * @brief Extents of component 'wrk_fou', leading dimension first.
     * @return All zeros when the component is unallocated.
     */
    std::array<int64_t, 2> wrk_fou_shape() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_curve_c_getcomp_wrk_fou(as<fitpack_curve_c>(), &raw, extents);
        return { extents[0], extents[1] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'wrk_fou'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through wrk_fou()(i, j) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> wrk_fou() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_curve_c_getcomp_wrk_fou(as<fitpack_curve_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "wrk_fou", FX_LEN_NAME - 1);
        descr.base_address = static_cast<void*>(raw);
        descr.type = getCFITypeFlag<double>();
        descr.elem_bytes = static_cast<FX_SIZE>(sizeof(double));
        descr.rank = static_cast<FX_RANK>(2);
        descr.is_pointer = true;   // non-owning: the object still owns the storage
        descr.is_slice = false;
        descr.attribute = static_cast<FX_ATTR>(FX_ATTR_POINTER);
        FX_SIZE stride = descr.elem_bytes;
        for (int k = 0; k < 2; ++k) {
            descr.dim[k].lower_bound = 0;   // C 0-based; Fortran lbound 1
            descr.dim[k].extent = static_cast<FX_SIZE>(extents[k]);
            descr.dim[k].stride_bytes = stride;
            stride *= descr.dim[k].extent;
        }
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    /**
     * @brief Deep copy of component 't' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> t_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_curve_c_getcomp_t(as<fitpack_curve_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 't'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through t()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> t() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_curve_c_getcomp_t(as<fitpack_curve_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "t", FX_LEN_NAME - 1);
        descr.base_address = static_cast<void*>(raw);
        descr.type = getCFITypeFlag<double>();
        descr.elem_bytes = static_cast<FX_SIZE>(sizeof(double));
        descr.rank = static_cast<FX_RANK>(1);
        descr.is_pointer = true;   // non-owning: the object still owns the storage
        descr.is_slice = false;
        descr.attribute = static_cast<FX_ATTR>(FX_ATTR_POINTER);
        FX_SIZE stride = descr.elem_bytes;
        for (int k = 0; k < 1; ++k) {
            descr.dim[k].lower_bound = 0;   // C 0-based; Fortran lbound 1
            descr.dim[k].extent = static_cast<FX_SIZE>(extents[k]);
            descr.dim[k].stride_bytes = stride;
            stride *= descr.dim[k].extent;
        }
        return fxArray<double>(descr);
    }
#endif // HAVE_FXARRAY

    // ===========================================================================================
    // Scalar Property Accessors
    // ===========================================================================================

    int32_t& m() {
        return *fitpack_curve_c_ref_m(as<fitpack_curve_c>());
    }
    const int32_t& m() const {
        return *fitpack_curve_c_ref_m(as<fitpack_curve_c>());
    }

    int32_t& order() {
        return *fitpack_curve_c_ref_order(as<fitpack_curve_c>());
    }
    const int32_t& order() const {
        return *fitpack_curve_c_ref_order(as<fitpack_curve_c>());
    }

    double& xleft() {
        return *fitpack_curve_c_ref_xleft(as<fitpack_curve_c>());
    }
    const double& xleft() const {
        return *fitpack_curve_c_ref_xleft(as<fitpack_curve_c>());
    }

    double& xright() {
        return *fitpack_curve_c_ref_xright(as<fitpack_curve_c>());
    }
    const double& xright() const {
        return *fitpack_curve_c_ref_xright(as<fitpack_curve_c>());
    }

    int32_t& nest() {
        return *fitpack_curve_c_ref_nest(as<fitpack_curve_c>());
    }
    const int32_t& nest() const {
        return *fitpack_curve_c_ref_nest(as<fitpack_curve_c>());
    }

    int32_t& bc() {
        return *fitpack_curve_c_ref_bc(as<fitpack_curve_c>());
    }
    const int32_t& bc() const {
        return *fitpack_curve_c_ref_bc(as<fitpack_curve_c>());
    }

    int32_t& knots() {
        return *fitpack_curve_c_ref_knots(as<fitpack_curve_c>());
    }
    const int32_t& knots() const {
        return *fitpack_curve_c_ref_knots(as<fitpack_curve_c>());
    }

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
    FP_REAL ddx(FP_REAL x, FP_SIZE order, FP_FLAG* ierr = nullptr) const { return dfdx(x, order, ierr); }

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

protected:
    explicit fpCurve(NoAlloc tag) : fxFitpackFitter(tag) {}

};

#endif /* FPCURVE_HPP_INCLUDED */
