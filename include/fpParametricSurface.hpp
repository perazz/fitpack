/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fpParametricSurface.hpp                                                                     **
 *   ** @brief Standalone C++ wrapper for fitpack_parametric_surface (no fortran-arrays dependency)
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FPPARAMETRICSURFACE_HPP_INCLUDED
#define FPPARAMETRICSURFACE_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_parametric_surface_c.h"
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

static_assert(sizeof(fitpack_parametric_surface_c) == sizeof(fitpack_fitter_c),
    "C descriptor layout mismatch: fitpack_parametric_surface_c vs fitpack_fitter_c");

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_parametric_surface
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_parametric_surface, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 * Extends fxFitpackFitter (Fortran: extends(fitpack_fitter))
 */
class fpParametricSurface : public fxFitpackFitter {
public:
    // ===========================================================================================
    // Constructors and Destructor
    // ===========================================================================================

    /**
     * @brief Default constructor - allocates new fitpack_parametric_surface
     */
    fpParametricSurface() : fxFitpackFitter(NoAlloc{}) {
        fitpack_parametric_surface_c_allocate(as<fitpack_parametric_surface_c>(), nullptr);
    }

    /**
     * @brief Destructor - deallocates if owned
     */
    ~fpParametricSurface() override {
        fitpack_parametric_surface_c_destroy(as<fitpack_parametric_surface_c>(), nullptr);
    }

    /**
     * @brief Copy constructor - deep copy
     */
    fpParametricSurface(const fpParametricSurface& other) : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_parametric_surface_c>() = fitpack_parametric_surface_c_null;
        fitpack_parametric_surface_c_copy(as<fitpack_parametric_surface_c>(), other.as<fitpack_parametric_surface_c>(), false, nullptr);
    }

    /**
     * @brief Copy assignment - deep copy
     */
    fpParametricSurface& operator=(const fpParametricSurface& other) {
        if (this != &other) {
            fitpack_parametric_surface_c_destroy(as<fitpack_parametric_surface_c>(), nullptr);
            fitpack_parametric_surface_c_copy(as<fitpack_parametric_surface_c>(), other.as<fitpack_parametric_surface_c>(), false, nullptr);
        }
        return *this;
    }

    /**
     * @brief Move constructor - transfer ownership
     */
    fpParametricSurface(fpParametricSurface&& other) noexcept : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_parametric_surface_c>() = fitpack_parametric_surface_c_null;
        fitpack_parametric_surface_c_move_alloc(as<fitpack_parametric_surface_c>(), other.as<fitpack_parametric_surface_c>(), nullptr);
    }

    /**
     * @brief Move assignment - transfer ownership
     */
    fpParametricSurface& operator=(fpParametricSurface&& other) noexcept {
        if (this != &other) {
            fitpack_parametric_surface_c_destroy(as<fitpack_parametric_surface_c>(), nullptr);
            fitpack_parametric_surface_c_move_alloc(as<fitpack_parametric_surface_c>(), other.as<fitpack_parametric_surface_c>(), nullptr);
        }
        return *this;
    }

    /**
     * @brief Construct from existing C wrapper (takes ownership if move=true)
     */
    explicit fpParametricSurface(fitpack_parametric_surface_c& c_wrapper, bool move = false) : fxFitpackFitter(NoAlloc{}) {
        if (move) {
            fitpack_parametric_surface_c_move_alloc(as<fitpack_parametric_surface_c>(), &c_wrapper, nullptr);
        } else {
            *as<fitpack_parametric_surface_c>() = fitpack_parametric_surface_c_null;
            fitpack_parametric_surface_c_copy(as<fitpack_parametric_surface_c>(), &c_wrapper, false, nullptr);
        }
    }

    /**
     * @brief Non-owning view ctor: bit-copy the C handle, mark
     * is_pointer=true. Used by parent classes' polymorphic accessors
     * to wrap a base-class handle as the correct concrete subtype.
     */
    explicit fpParametricSurface(fitpack_parametric_surface_c& c_wrapper, ViewTag) : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_parametric_surface_c>() = c_wrapper;
        as<fitpack_parametric_surface_c>()->is_pointer = true;
    }

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    bool is_allocated() const override {
        return as<fitpack_parametric_surface_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    bool is_pointer() const override {
        return as<fitpack_parametric_surface_c>()->is_pointer;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_parametric_surface_c").
     */
    const char* c_type_name() const override {
        return fitpack_parametric_surface_c_c_type_name(*as<fitpack_parametric_surface_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    const char* cpp_type_name() const override {
        return "fpParametricSurface";
    }

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_parametric_surface_c& c_handle() {
        return *as<fitpack_parametric_surface_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_parametric_surface_c& c_handle() const {
        return *as<fitpack_parametric_surface_c>();
    }

    /**
     * @brief Construct an empty (non-owning) wrapper, for use as a view target.
     *
     * The returned object has a null handle (cptr==nullptr, is_pointer==false)
     * and does NOT allocate a Fortran object. Intended as a fill target for
     * member-view accessors emitted on parent types — the parent populates
     * the handle with c_loc() into its own storage and sets is_pointer=true.
     */
    static fpParametricSurface make_view() {
        return fpParametricSurface(NoAlloc{});
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
     * @brief fit
     */
    virtual int32_t fit(int32_t n, double* smoothing = nullptr, bool* periodic = nullptr, bool* keep_knots = nullptr) const {
        return fitpack_parametric_surface_c_fit(as<fitpack_parametric_surface_c>(), smoothing, n, periodic, keep_knots);
    }


    virtual int32_t interpolate(bool* reset_knots = nullptr) const {
        return fitpack_parametric_surface_c_interpolate(as<fitpack_parametric_surface_c>(), reset_knots);
    }

    /**
     * @brief least_squares
     */
    virtual int32_t least_squares(int32_t n, double* u_knots = nullptr, double* v_knots = nullptr, double* smoothing = nullptr, bool* reset_knots = nullptr) const {
        return fitpack_parametric_surface_c_least_squares(as<fitpack_parametric_surface_c>(), n, u_knots, v_knots, smoothing, reset_knots);
    }


    /**
     * @brief eval
     */
    virtual std::vector<double> eval(double u, double v, int32_t n_result, int32_t* ierr = nullptr) {
        std::vector<double> result(n_result);
        fitpack_parametric_surface_c_surf_eval_one(as<fitpack_parametric_surface_c>(), u, v, ierr, result.data(), n_result);
        return result;
    }


    /**
     * @brief eval
     */
    virtual std::vector<double> eval(std::vector<double>& u, std::vector<double>& v, int32_t n_result, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(u.size());
        std::vector<double> result(n_result);
        fitpack_parametric_surface_c_surf_eval_grid(as<fitpack_parametric_surface_c>(), n, u.data(), v.data(), ierr, result.data(), n_result);
        return result;
    }


    int32_t comm_size() const override {
        return fitpack_parametric_surface_c_comm_size(as<fitpack_parametric_surface_c>());
    }

    /**
     * @brief comm_pack
     */
    void comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_parametric_surface_c_comm_pack(as<fitpack_parametric_surface_c>(), n, buffer.data());
    }


    /**
     * @brief comm_expand
     */
    void comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_parametric_surface_c_comm_expand(as<fitpack_parametric_surface_c>(), n, buffer.data());
    }


    double mse() const override {
        return fitpack_parametric_surface_c_mse(as<fitpack_parametric_surface_c>());
    }

    int32_t core_comm_size() const override {
        return fitpack_parametric_surface_c_core_comm_size(as<fitpack_parametric_surface_c>());
    }

    /**
     * @brief core_comm_pack
     */
    void core_comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_parametric_surface_c_core_comm_pack(as<fitpack_parametric_surface_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    void core_comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_parametric_surface_c_core_comm_expand(as<fitpack_parametric_surface_c>(), n, buffer.data());
    }


    void destroy_base() override {
        fitpack_parametric_surface_c_destroy_base(as<fitpack_parametric_surface_c>());
    }

    // ===========================================================================================
    // Component Array Accessors
    // ===========================================================================================

    /**
     * @brief Deep copy of component 'u' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> u_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_parametric_surface_c_getcomp_u(as<fitpack_parametric_surface_c>(), &raw, extents);
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
        fitpack_parametric_surface_c_getcomp_u(as<fitpack_parametric_surface_c>(), &raw, extents);
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
     * @brief Deep copy of component 'v' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> v_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_parametric_surface_c_getcomp_v(as<fitpack_parametric_surface_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'v'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through v()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> v() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_parametric_surface_c_getcomp_v(as<fitpack_parametric_surface_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "v", FX_LEN_NAME - 1);
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
     * @brief Deep copy of component 'z' as a std::vector.
     *
     * The rank-3 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with z_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> z_vector() const {
        double* raw = nullptr;
        int64_t extents[3] = {0, 0, 0};
        fitpack_parametric_surface_c_getcomp_z(as<fitpack_parametric_surface_c>(), &raw, extents);
        const int64_t total = extents[0] * extents[1] * extents[2];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

    /**
     * @brief Extents of component 'z', leading dimension first.
     * @return All zeros when the component is unallocated.
     */
    std::array<int64_t, 3> z_shape() const {
        double* raw = nullptr;
        int64_t extents[3] = {0, 0, 0};
        fitpack_parametric_surface_c_getcomp_z(as<fitpack_parametric_surface_c>(), &raw, extents);
        return { extents[0], extents[1], extents[2] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'z'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through z()(i, j) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> z() const {
        double* raw = nullptr;
        int64_t extents[3] = {0, 0, 0};
        fitpack_parametric_surface_c_getcomp_z(as<fitpack_parametric_surface_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "z", FX_LEN_NAME - 1);
        descr.base_address = static_cast<void*>(raw);
        descr.type = getCFITypeFlag<double>();
        descr.elem_bytes = static_cast<FX_SIZE>(sizeof(double));
        descr.rank = static_cast<FX_RANK>(3);
        descr.is_pointer = true;   // non-owning: the object still owns the storage
        descr.is_slice = false;
        descr.attribute = static_cast<FX_ATTR>(FX_ATTR_POINTER);
        FX_SIZE stride = descr.elem_bytes;
        for (int k = 0; k < 3; ++k) {
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
     *
     * The rank-2 component is flattened COLUMN-MAJOR (Fortran order):
     * element (i, j, ...) — zero-based — is at index i + j * extents[0] + ...
     * Query the extents with t_shape().
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> t_vector() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_parametric_surface_c_getcomp_t(as<fitpack_parametric_surface_c>(), &raw, extents);
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
        fitpack_parametric_surface_c_getcomp_t(as<fitpack_parametric_surface_c>(), &raw, extents);
        return { extents[0], extents[1] };
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 't'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through t()(i, j) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> t() const {
        double* raw = nullptr;
        int64_t extents[2] = {0, 0};
        fitpack_parametric_surface_c_getcomp_t(as<fitpack_parametric_surface_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "t", FX_LEN_NAME - 1);
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

    // ===========================================================================================
    // Scalar Property Accessors
    // ===========================================================================================

    int32_t& idim() {
        return *fitpack_parametric_surface_c_ref_idim(as<fitpack_parametric_surface_c>());
    }
    const int32_t& idim() const {
        return *fitpack_parametric_surface_c_ref_idim(as<fitpack_parametric_surface_c>());
    }

    int32_t& nmax() {
        return *fitpack_parametric_surface_c_ref_nmax(as<fitpack_parametric_surface_c>());
    }
    const int32_t& nmax() const {
        return *fitpack_parametric_surface_c_ref_nmax(as<fitpack_parametric_surface_c>());
    }

protected:
    explicit fpParametricSurface(NoAlloc tag) : fxFitpackFitter(tag) {}

};

#endif /* FPPARAMETRICSURFACE_HPP_INCLUDED */
