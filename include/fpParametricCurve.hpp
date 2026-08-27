/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fpParametricCurve.hpp                                                                     **
 *   ** @brief Standalone C++ wrapper for fitpack_parametric_curve (no fortran-arrays dependency)
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FPPARAMETRICCURVE_HPP_INCLUDED
#define FPPARAMETRICCURVE_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_parametric_curve_c.h"
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

static_assert(sizeof(fitpack_parametric_curve_c) == sizeof(fitpack_fitter_c),
    "C descriptor layout mismatch: fitpack_parametric_curve_c vs fitpack_fitter_c");

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_parametric_curve
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_parametric_curve, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 * Extends fxFitpackFitter (Fortran: extends(fitpack_fitter))
 */
class fpParametricCurve : public fxFitpackFitter {
public:
    // ===========================================================================================
    // Constructors and Destructor
    // ===========================================================================================

    /**
     * @brief Default constructor - allocates new fitpack_parametric_curve
     */
    fpParametricCurve() : fxFitpackFitter(NoAlloc{}) {
        fitpack_parametric_curve_c_allocate(as<fitpack_parametric_curve_c>(), nullptr);
    }

    /**
     * @brief Destructor - deallocates if owned
     */
    ~fpParametricCurve() override {
        fitpack_parametric_curve_c_destroy(as<fitpack_parametric_curve_c>(), nullptr);
    }

    /**
     * @brief Copy constructor - deep copy
     */
    fpParametricCurve(const fpParametricCurve& other) : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_parametric_curve_c>() = fitpack_parametric_curve_c_null;
        fitpack_parametric_curve_c_copy(as<fitpack_parametric_curve_c>(), other.as<fitpack_parametric_curve_c>(), false, nullptr);
    }

    /**
     * @brief Copy assignment - deep copy
     */
    fpParametricCurve& operator=(const fpParametricCurve& other) {
        if (this != &other) {
            fitpack_parametric_curve_c_destroy(as<fitpack_parametric_curve_c>(), nullptr);
            fitpack_parametric_curve_c_copy(as<fitpack_parametric_curve_c>(), other.as<fitpack_parametric_curve_c>(), false, nullptr);
        }
        return *this;
    }

    /**
     * @brief Move constructor - transfer ownership
     */
    fpParametricCurve(fpParametricCurve&& other) noexcept : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_parametric_curve_c>() = fitpack_parametric_curve_c_null;
        fitpack_parametric_curve_c_move_alloc(as<fitpack_parametric_curve_c>(), other.as<fitpack_parametric_curve_c>(), nullptr);
    }

    /**
     * @brief Move assignment - transfer ownership
     */
    fpParametricCurve& operator=(fpParametricCurve&& other) noexcept {
        if (this != &other) {
            fitpack_parametric_curve_c_destroy(as<fitpack_parametric_curve_c>(), nullptr);
            fitpack_parametric_curve_c_move_alloc(as<fitpack_parametric_curve_c>(), other.as<fitpack_parametric_curve_c>(), nullptr);
        }
        return *this;
    }

    /**
     * @brief Construct from existing C wrapper (takes ownership if move=true)
     */
    explicit fpParametricCurve(fitpack_parametric_curve_c& c_wrapper, bool move = false) : fxFitpackFitter(NoAlloc{}) {
        if (move) {
            fitpack_parametric_curve_c_move_alloc(as<fitpack_parametric_curve_c>(), &c_wrapper, nullptr);
        } else {
            *as<fitpack_parametric_curve_c>() = fitpack_parametric_curve_c_null;
            fitpack_parametric_curve_c_copy(as<fitpack_parametric_curve_c>(), &c_wrapper, false, nullptr);
        }
    }

    /**
     * @brief Non-owning view ctor: bit-copy the C handle, mark
     * is_pointer=true. Used by parent classes' polymorphic accessors
     * to wrap a base-class handle as the correct concrete subtype.
     */
    explicit fpParametricCurve(fitpack_parametric_curve_c& c_wrapper, ViewTag) : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_parametric_curve_c>() = c_wrapper;
        as<fitpack_parametric_curve_c>()->is_pointer = true;
    }

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    bool is_allocated() const override {
        return as<fitpack_parametric_curve_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    bool is_pointer() const override {
        return as<fitpack_parametric_curve_c>()->is_pointer;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_parametric_curve_c").
     */
    const char* c_type_name() const override {
        return fitpack_parametric_curve_c_c_type_name(*as<fitpack_parametric_curve_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    const char* cpp_type_name() const override {
        return "fpParametricCurve";
    }

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_parametric_curve_c& c_handle() {
        return *as<fitpack_parametric_curve_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_parametric_curve_c& c_handle() const {
        return *as<fitpack_parametric_curve_c>();
    }

    /**
     * @brief Construct an empty (non-owning) wrapper, for use as a view target.
     *
     * The returned object has a null handle (cptr==nullptr, is_pointer==false)
     * and does NOT allocate a Fortran object. Intended as a fill target for
     * member-view accessors emitted on parent types — the parent populates
     * the handle with c_loc() into its own storage and sets is_pointer=true.
     */
    static fpParametricCurve make_view() {
        return fpParametricCurve(NoAlloc{});
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
    virtual void new_points(int32_t x_n1, int32_t x_n2, std::vector<double>& x, double* u = nullptr, double* w = nullptr) {
        fitpack_parametric_curve_c_new_points(as<fitpack_parametric_curve_c>(), x_n1, x_n2, x.data(), u, w);
    }


    virtual void set_default_parameters() {
        fitpack_parametric_curve_c_set_default_parameters(as<fitpack_parametric_curve_c>());
    }

    /**
     * @brief new_fit
     */
    virtual int32_t new_fit(int32_t x_n1, int32_t x_n2, std::vector<double>& x, double* u = nullptr, double* w = nullptr, double* smoothing = nullptr, int32_t* order = nullptr) const {
        return fitpack_parametric_curve_c_new_fit(as<fitpack_parametric_curve_c>(), x_n1, x_n2, x.data(), u, w, smoothing, order);
    }


    virtual int32_t fit(double* smoothing = nullptr, int32_t* order = nullptr, bool* keep_knots = nullptr) const {
        return fitpack_parametric_curve_c_fit(as<fitpack_parametric_curve_c>(), smoothing, order, keep_knots);
    }

    virtual int32_t interpolate(int32_t* order = nullptr, bool* reset_knots = nullptr) const {
        return fitpack_parametric_curve_c_interpolate(as<fitpack_parametric_curve_c>(), order, reset_knots);
    }

    virtual int32_t least_squares(double* smoothing = nullptr, bool* reset_knots = nullptr) const {
        return fitpack_parametric_curve_c_least_squares(as<fitpack_parametric_curve_c>(), smoothing, reset_knots);
    }

    /**
     * @brief eval_one
     */
    virtual std::vector<double> eval_one(double u, int32_t n_result, int32_t* ierr = nullptr) {
        std::vector<double> result(n_result);
        fitpack_parametric_curve_c_eval_one(as<fitpack_parametric_curve_c>(), u, ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief eval_many
     */
    virtual std::vector<double> eval_many(std::vector<double>& u, int32_t n_result, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(u.size());
        std::vector<double> result(n_result);
        fitpack_parametric_curve_c_eval_many(as<fitpack_parametric_curve_c>(), n, u.data(), ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx
     */
    virtual std::vector<double> dfdx(double u, int32_t order, int32_t n_result, int32_t* ierr = nullptr) {
        std::vector<double> result(n_result);
        fitpack_parametric_curve_c_curve_derivative(as<fitpack_parametric_curve_c>(), u, order, ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx
     */
    virtual std::vector<double> dfdx(std::vector<double>& u, int32_t order, int32_t n_result, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(u.size());
        std::vector<double> result(n_result);
        fitpack_parametric_curve_c_curve_derivatives(as<fitpack_parametric_curve_c>(), n, u.data(), order, ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief dfdx_all
     */
    virtual std::vector<double> dfdx_all(double u, int32_t n_result, int32_t* ierr = nullptr) {
        std::vector<double> result(n_result);
        fitpack_parametric_curve_c_curve_all_derivatives(as<fitpack_parametric_curve_c>(), u, ierr, result.data(), n_result);
        return result;
    }


    int32_t comm_size() const override {
        return fitpack_parametric_curve_c_comm_size(as<fitpack_parametric_curve_c>());
    }

    /**
     * @brief comm_pack
     */
    void comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_parametric_curve_c_comm_pack(as<fitpack_parametric_curve_c>(), n, buffer.data());
    }


    /**
     * @brief comm_expand
     */
    void comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_parametric_curve_c_comm_expand(as<fitpack_parametric_curve_c>(), n, buffer.data());
    }


    double mse() const override {
        return fitpack_parametric_curve_c_mse(as<fitpack_parametric_curve_c>());
    }

    int32_t core_comm_size() const override {
        return fitpack_parametric_curve_c_core_comm_size(as<fitpack_parametric_curve_c>());
    }

    /**
     * @brief core_comm_pack
     */
    void core_comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_parametric_curve_c_core_comm_pack(as<fitpack_parametric_curve_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    void core_comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_parametric_curve_c_core_comm_expand(as<fitpack_parametric_curve_c>(), n, buffer.data());
    }


    void destroy_base() override {
        fitpack_parametric_curve_c_destroy_base(as<fitpack_parametric_curve_c>());
    }

    // ===========================================================================================
    // Component Array Accessors
    // ===========================================================================================

    /**
     * @brief Deep copy of component 'x' as a std::vector.
     *
     * The rank-2 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with x_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> x_vector() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_parametric_curve_c_getcomp_x(as<fitpack_parametric_curve_c>(), &raw, extents);
        const int64_t total = extents[0] * extents[1];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

    /**
     * @brief Extents of component 'x', leading dimension first.
     * @return All zeros when the component is unallocated.
     */
    std::array<int64_t, 2> x_shape() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_parametric_curve_c_getcomp_x(as<fitpack_parametric_curve_c>(), &raw, extents);
        return { extents[0], extents[1] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'x'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through x()(i, j) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> x() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_parametric_curve_c_getcomp_x(as<fitpack_parametric_curve_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "x", FX_LEN_NAME - 1);
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
     * @brief Deep copy of component 'u' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> u_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_parametric_curve_c_getcomp_u(as<fitpack_parametric_curve_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'u'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through u()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> u() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_parametric_curve_c_getcomp_u(as<fitpack_parametric_curve_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "u", FX_LEN_NAME - 1);
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
        fitpack_parametric_curve_c_getcomp_sp(as<fitpack_parametric_curve_c>(), &raw, extents);
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
        fitpack_parametric_curve_c_getcomp_sp(as<fitpack_parametric_curve_c>(), &raw, extents);
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
        fitpack_parametric_curve_c_getcomp_w(as<fitpack_parametric_curve_c>(), &raw, extents);
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
        fitpack_parametric_curve_c_getcomp_w(as<fitpack_parametric_curve_c>(), &raw, extents);
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
     * @brief Deep copy of component 'dd' as a std::vector.
     *
     * The rank-2 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with dd_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> dd_vector() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_parametric_curve_c_getcomp_dd(as<fitpack_parametric_curve_c>(), &raw, extents);
        const int64_t total = extents[0] * extents[1];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

    /**
     * @brief Extents of component 'dd', leading dimension first.
     * @return All zeros when the component is unallocated.
     */
    std::array<int64_t, 2> dd_shape() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_parametric_curve_c_getcomp_dd(as<fitpack_parametric_curve_c>(), &raw, extents);
        return { extents[0], extents[1] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'dd'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through dd()(i, j) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> dd() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_parametric_curve_c_getcomp_dd(as<fitpack_parametric_curve_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "dd", FX_LEN_NAME - 1);
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
        fitpack_parametric_curve_c_getcomp_t(as<fitpack_parametric_curve_c>(), &raw, extents);
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
        fitpack_parametric_curve_c_getcomp_t(as<fitpack_parametric_curve_c>(), &raw, extents);
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
        return *fitpack_parametric_curve_c_ref_m(as<fitpack_parametric_curve_c>());
    }
    const int32_t& m() const {
        return *fitpack_parametric_curve_c_ref_m(as<fitpack_parametric_curve_c>());
    }

    int32_t& idim() {
        return *fitpack_parametric_curve_c_ref_idim(as<fitpack_parametric_curve_c>());
    }
    const int32_t& idim() const {
        return *fitpack_parametric_curve_c_ref_idim(as<fitpack_parametric_curve_c>());
    }

    bool has_params() const {
        return fitpack_parametric_curve_c_get_has_params(as<fitpack_parametric_curve_c>());
    }
    void set_has_params(bool value) {
        fitpack_parametric_curve_c_set_has_params(as<fitpack_parametric_curve_c>(), value);
    }

    int32_t& order() {
        return *fitpack_parametric_curve_c_ref_order(as<fitpack_parametric_curve_c>());
    }
    const int32_t& order() const {
        return *fitpack_parametric_curve_c_ref_order(as<fitpack_parametric_curve_c>());
    }

    double& ubegin() {
        return *fitpack_parametric_curve_c_ref_ubegin(as<fitpack_parametric_curve_c>());
    }
    const double& ubegin() const {
        return *fitpack_parametric_curve_c_ref_ubegin(as<fitpack_parametric_curve_c>());
    }

    double& uend() {
        return *fitpack_parametric_curve_c_ref_uend(as<fitpack_parametric_curve_c>());
    }
    const double& uend() const {
        return *fitpack_parametric_curve_c_ref_uend(as<fitpack_parametric_curve_c>());
    }

    int32_t& nest() {
        return *fitpack_parametric_curve_c_ref_nest(as<fitpack_parametric_curve_c>());
    }
    const int32_t& nest() const {
        return *fitpack_parametric_curve_c_ref_nest(as<fitpack_parametric_curve_c>());
    }

    int32_t& knots() {
        return *fitpack_parametric_curve_c_ref_knots(as<fitpack_parametric_curve_c>());
    }
    const int32_t& knots() const {
        return *fitpack_parametric_curve_c_ref_knots(as<fitpack_parametric_curve_c>());
    }

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

protected:
    explicit fpParametricCurve(NoAlloc tag) : fxFitpackFitter(tag) {}

};

#endif /* FPPARAMETRICCURVE_HPP_INCLUDED */
