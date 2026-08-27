/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fpSphere.hpp                                                                     **
 *   ** @brief Standalone C++ wrapper for fitpack_sphere (no fortran-arrays dependency)
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FPSPHERE_HPP_INCLUDED
#define FPSPHERE_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_sphere_c.h"
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

static_assert(sizeof(fitpack_sphere_c) == sizeof(fitpack_fitter_c),
    "C descriptor layout mismatch: fitpack_sphere_c vs fitpack_fitter_c");

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_sphere
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_sphere, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 * Extends fxFitpackFitter (Fortran: extends(fitpack_fitter))
 */
class fpSphere : public fxFitpackFitter {
public:
    // ===========================================================================================
    // Constructors and Destructor
    // ===========================================================================================

    /**
     * @brief Default constructor - allocates new fitpack_sphere
     */
    fpSphere() : fxFitpackFitter(NoAlloc{}) {
        fitpack_sphere_c_allocate(as<fitpack_sphere_c>(), nullptr);
    }

    /**
     * @brief Destructor - deallocates if owned
     */
    ~fpSphere() override {
        fitpack_sphere_c_destroy(as<fitpack_sphere_c>(), nullptr);
    }

    /**
     * @brief Copy constructor - deep copy
     */
    fpSphere(const fpSphere& other) : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_sphere_c>() = fitpack_sphere_c_null;
        fitpack_sphere_c_copy(as<fitpack_sphere_c>(), other.as<fitpack_sphere_c>(), false, nullptr);
    }

    /**
     * @brief Copy assignment - deep copy
     */
    fpSphere& operator=(const fpSphere& other) {
        if (this != &other) {
            fitpack_sphere_c_destroy(as<fitpack_sphere_c>(), nullptr);
            fitpack_sphere_c_copy(as<fitpack_sphere_c>(), other.as<fitpack_sphere_c>(), false, nullptr);
        }
        return *this;
    }

    /**
     * @brief Move constructor - transfer ownership
     */
    fpSphere(fpSphere&& other) noexcept : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_sphere_c>() = fitpack_sphere_c_null;
        fitpack_sphere_c_move_alloc(as<fitpack_sphere_c>(), other.as<fitpack_sphere_c>(), nullptr);
    }

    /**
     * @brief Move assignment - transfer ownership
     */
    fpSphere& operator=(fpSphere&& other) noexcept {
        if (this != &other) {
            fitpack_sphere_c_destroy(as<fitpack_sphere_c>(), nullptr);
            fitpack_sphere_c_move_alloc(as<fitpack_sphere_c>(), other.as<fitpack_sphere_c>(), nullptr);
        }
        return *this;
    }

    /**
     * @brief Construct from existing C wrapper (takes ownership if move=true)
     */
    explicit fpSphere(fitpack_sphere_c& c_wrapper, bool move = false) : fxFitpackFitter(NoAlloc{}) {
        if (move) {
            fitpack_sphere_c_move_alloc(as<fitpack_sphere_c>(), &c_wrapper, nullptr);
        } else {
            *as<fitpack_sphere_c>() = fitpack_sphere_c_null;
            fitpack_sphere_c_copy(as<fitpack_sphere_c>(), &c_wrapper, false, nullptr);
        }
    }

    /**
     * @brief Non-owning view ctor: bit-copy the C handle, mark
     * is_pointer=true. Used by parent classes' polymorphic accessors
     * to wrap a base-class handle as the correct concrete subtype.
     */
    explicit fpSphere(fitpack_sphere_c& c_wrapper, ViewTag) : fxFitpackFitter(NoAlloc{}) {
        *as<fitpack_sphere_c>() = c_wrapper;
        as<fitpack_sphere_c>()->is_pointer = true;
    }

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    bool is_allocated() const override {
        return as<fitpack_sphere_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    bool is_pointer() const override {
        return as<fitpack_sphere_c>()->is_pointer;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_sphere_c").
     */
    const char* c_type_name() const override {
        return fitpack_sphere_c_c_type_name(*as<fitpack_sphere_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    const char* cpp_type_name() const override {
        return "fpSphere";
    }

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_sphere_c& c_handle() {
        return *as<fitpack_sphere_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_sphere_c& c_handle() const {
        return *as<fitpack_sphere_c>();
    }

    /**
     * @brief Construct an empty (non-owning) wrapper, for use as a view target.
     *
     * The returned object has a null handle (cptr==nullptr, is_pointer==false)
     * and does NOT allocate a Fortran object. Intended as a fill target for
     * member-view accessors emitted on parent types — the parent populates
     * the handle with c_loc() into its own storage and sets is_pointer=true.
     */
    static fpSphere make_view() {
        return fpSphere(NoAlloc{});
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
    virtual void new_points(std::vector<double>& theta, std::vector<double>& phi, std::vector<double>& r, double* w = nullptr) {
        int32_t n = static_cast<int32_t>(theta.size());
        fitpack_sphere_c_new_points(as<fitpack_sphere_c>(), n, theta.data(), phi.data(), r.data(), w);
    }


    /**
     * @brief new_fit
     */
    virtual int32_t new_fit(std::vector<double>& theta, std::vector<double>& phi, std::vector<double>& r, double* w = nullptr, double* smoothing = nullptr) const {
        int32_t n = static_cast<int32_t>(theta.size());
        return fitpack_sphere_c_new_fit(as<fitpack_sphere_c>(), n, theta.data(), phi.data(), r.data(), w, smoothing);
    }


    virtual int32_t fit(double* smoothing = nullptr, bool* keep_knots = nullptr) const {
        return fitpack_sphere_c_fit(as<fitpack_sphere_c>(), smoothing, keep_knots);
    }

    virtual int32_t least_squares(double* smoothing = nullptr, bool* reset_knots = nullptr) const {
        return fitpack_sphere_c_least_squares(as<fitpack_sphere_c>(), smoothing, reset_knots);
    }

    virtual int32_t interpolate(bool* reset_knots = nullptr) const {
        return fitpack_sphere_c_interpolate(as<fitpack_sphere_c>(), reset_knots);
    }

    virtual double eval(double theta, double phi, int32_t* ierr = nullptr) const {
        return fitpack_sphere_c_sphere_eval_one(as<fitpack_sphere_c>(), theta, phi, ierr);
    }

    /**
     * @brief eval
     */
    virtual std::vector<double> eval(std::vector<double>& theta, std::vector<double>& phi, int32_t* ierr = nullptr) {
        int32_t n = static_cast<int32_t>(theta.size());
        std::vector<double> result(n * n);
        int32_t n_result = 0;
        fitpack_sphere_c_sphere_eval_many(as<fitpack_sphere_c>(), n, theta.data(), phi.data(), ierr, result.data(), &n_result, n * n);
        result.resize(n_result);
        return result;
    }


    int32_t comm_size() const override {
        return fitpack_sphere_c_comm_size(as<fitpack_sphere_c>());
    }

    /**
     * @brief comm_pack
     */
    void comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_sphere_c_comm_pack(as<fitpack_sphere_c>(), n, buffer.data());
    }


    /**
     * @brief comm_expand
     */
    void comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_sphere_c_comm_expand(as<fitpack_sphere_c>(), n, buffer.data());
    }


    double mse() const override {
        return fitpack_sphere_c_mse(as<fitpack_sphere_c>());
    }

    int32_t core_comm_size() const override {
        return fitpack_sphere_c_core_comm_size(as<fitpack_sphere_c>());
    }

    /**
     * @brief core_comm_pack
     */
    void core_comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_sphere_c_core_comm_pack(as<fitpack_sphere_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    void core_comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_sphere_c_core_comm_expand(as<fitpack_sphere_c>(), n, buffer.data());
    }


    void destroy_base() override {
        fitpack_sphere_c_destroy_base(as<fitpack_sphere_c>());
    }

    // ===========================================================================================
    // Component Array Accessors
    // ===========================================================================================

    /**
     * @brief Deep copy of component 'theta' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> theta_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_sphere_c_getcomp_theta(as<fitpack_sphere_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'theta'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through theta()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> theta() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_sphere_c_getcomp_theta(as<fitpack_sphere_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "theta", FX_LEN_NAME - 1);
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
     * @brief Deep copy of component 'phi' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> phi_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_sphere_c_getcomp_phi(as<fitpack_sphere_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'phi'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through phi()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> phi() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_sphere_c_getcomp_phi(as<fitpack_sphere_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "phi", FX_LEN_NAME - 1);
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
     * @brief Deep copy of component 'r' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> r_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_sphere_c_getcomp_r(as<fitpack_sphere_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'r'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through r()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> r() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_sphere_c_getcomp_r(as<fitpack_sphere_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "r", FX_LEN_NAME - 1);
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
        fitpack_sphere_c_getcomp_w(as<fitpack_sphere_c>(), &raw, extents);
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
        fitpack_sphere_c_getcomp_w(as<fitpack_sphere_c>(), &raw, extents);
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
     * @brief Deep copy of component 'wrk2' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> wrk2_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_sphere_c_getcomp_wrk2(as<fitpack_sphere_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'wrk2'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through wrk2()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> wrk2() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_sphere_c_getcomp_wrk2(as<fitpack_sphere_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "wrk2", FX_LEN_NAME - 1);
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
        fitpack_sphere_c_getcomp_t(as<fitpack_sphere_c>(), &raw, extents);
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
        fitpack_sphere_c_getcomp_t(as<fitpack_sphere_c>(), &raw, extents);
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
        fitpack_sphere_c_getcomp_t(as<fitpack_sphere_c>(), &raw, extents);
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

    int32_t& m() {
        return *fitpack_sphere_c_ref_m(as<fitpack_sphere_c>());
    }
    const int32_t& m() const {
        return *fitpack_sphere_c_ref_m(as<fitpack_sphere_c>());
    }

    int32_t& lwrk2() {
        return *fitpack_sphere_c_ref_lwrk2(as<fitpack_sphere_c>());
    }
    const int32_t& lwrk2() const {
        return *fitpack_sphere_c_ref_lwrk2(as<fitpack_sphere_c>());
    }

    int32_t& nmax() {
        return *fitpack_sphere_c_ref_nmax(as<fitpack_sphere_c>());
    }
    const int32_t& nmax() const {
        return *fitpack_sphere_c_ref_nmax(as<fitpack_sphere_c>());
    }

protected:
    explicit fpSphere(NoAlloc tag) : fxFitpackFitter(tag) {}

};

#endif /* FPSPHERE_HPP_INCLUDED */
