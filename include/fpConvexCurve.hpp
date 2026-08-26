/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fpConvexCurve.hpp                                                                     **
 *   ** @brief Standalone C++ wrapper for fitpack_convex_curve (no fortran-arrays dependency)
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FPCONVEXCURVE_HPP_INCLUDED
#define FPCONVEXCURVE_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_convex_curve_c.h"
#include "fpCurve.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <variant>
#include <optional>
#include <cstring>

static_assert(sizeof(fitpack_convex_curve_c) == sizeof(fitpack_curve_c),
    "C descriptor layout mismatch: fitpack_convex_curve_c vs fitpack_curve_c");

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_convex_curve
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_convex_curve, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 * Extends fpCurve (Fortran: extends(fitpack_curve))
 */
class fpConvexCurve : public fpCurve {
public:
    // ===========================================================================================
    // Constructors and Destructor
    // ===========================================================================================

    /**
     * @brief Default constructor - allocates new fitpack_convex_curve
     */
    fpConvexCurve() : fpCurve(NoAlloc{}) {
        fitpack_convex_curve_c_allocate(as<fitpack_convex_curve_c>(), nullptr);
    }

    /**
     * @brief Destructor - deallocates if owned
     */
    ~fpConvexCurve() override {
        fitpack_convex_curve_c_destroy(as<fitpack_convex_curve_c>(), nullptr);
    }

    /**
     * @brief Copy constructor - deep copy
     */
    fpConvexCurve(const fpConvexCurve& other) : fpCurve(NoAlloc{}) {
        *as<fitpack_convex_curve_c>() = fitpack_convex_curve_c_null;
        fitpack_convex_curve_c_copy(as<fitpack_convex_curve_c>(), other.as<fitpack_convex_curve_c>(), false, nullptr);
    }

    /**
     * @brief Copy assignment - deep copy
     */
    fpConvexCurve& operator=(const fpConvexCurve& other) {
        if (this != &other) {
            fitpack_convex_curve_c_destroy(as<fitpack_convex_curve_c>(), nullptr);
            fitpack_convex_curve_c_copy(as<fitpack_convex_curve_c>(), other.as<fitpack_convex_curve_c>(), false, nullptr);
        }
        return *this;
    }

    /**
     * @brief Move constructor - transfer ownership
     */
    fpConvexCurve(fpConvexCurve&& other) noexcept : fpCurve(NoAlloc{}) {
        *as<fitpack_convex_curve_c>() = fitpack_convex_curve_c_null;
        fitpack_convex_curve_c_move_alloc(as<fitpack_convex_curve_c>(), other.as<fitpack_convex_curve_c>(), nullptr);
    }

    /**
     * @brief Move assignment - transfer ownership
     */
    fpConvexCurve& operator=(fpConvexCurve&& other) noexcept {
        if (this != &other) {
            fitpack_convex_curve_c_destroy(as<fitpack_convex_curve_c>(), nullptr);
            fitpack_convex_curve_c_move_alloc(as<fitpack_convex_curve_c>(), other.as<fitpack_convex_curve_c>(), nullptr);
        }
        return *this;
    }

    /**
     * @brief Construct from existing C wrapper (takes ownership if move=true)
     */
    explicit fpConvexCurve(fitpack_convex_curve_c& c_wrapper, bool move = false) : fpCurve(NoAlloc{}) {
        if (move) {
            fitpack_convex_curve_c_move_alloc(as<fitpack_convex_curve_c>(), &c_wrapper, nullptr);
        } else {
            *as<fitpack_convex_curve_c>() = fitpack_convex_curve_c_null;
            fitpack_convex_curve_c_copy(as<fitpack_convex_curve_c>(), &c_wrapper, false, nullptr);
        }
    }

    /**
     * @brief Non-owning view ctor: bit-copy the C handle, mark
     * is_pointer=true. Used by parent classes' polymorphic accessors
     * to wrap a base-class handle as the correct concrete subtype.
     */
    explicit fpConvexCurve(fitpack_convex_curve_c& c_wrapper, ViewTag) : fpCurve(NoAlloc{}) {
        *as<fitpack_convex_curve_c>() = c_wrapper;
        as<fitpack_convex_curve_c>()->is_pointer = true;
    }

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    bool is_allocated() const override {
        return as<fitpack_convex_curve_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    bool is_pointer() const override {
        return as<fitpack_convex_curve_c>()->is_pointer;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_convex_curve_c").
     */
    const char* c_type_name() const override {
        return fitpack_convex_curve_c_c_type_name(*as<fitpack_convex_curve_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    const char* cpp_type_name() const override {
        return "fpConvexCurve";
    }

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_convex_curve_c& c_handle() {
        return *as<fitpack_convex_curve_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_convex_curve_c& c_handle() const {
        return *as<fitpack_convex_curve_c>();
    }

    /**
     * @brief Construct an empty (non-owning) wrapper, for use as a view target.
     *
     * The returned object has a null handle (cptr==nullptr, is_pointer==false)
     * and does NOT allocate a Fortran object. Intended as a fill target for
     * member-view accessors emitted on parent types — the parent populates
     * the handle with c_loc() into its own storage and sets is_pointer=true.
     */
    static fpConvexCurve make_view() {
        return fpConvexCurve(NoAlloc{});
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
        fitpack_convex_curve_c_new_points(as<fitpack_convex_curve_c>(), n, x.data(), y.data(), w);
    }


    int32_t fit(double* smoothing = nullptr, int32_t* order = nullptr, bool* keep_knots = nullptr) const override {
        return fitpack_convex_curve_c_fit(as<fitpack_convex_curve_c>(), smoothing, order, keep_knots);
    }

    int32_t least_squares(double* smoothing = nullptr, bool* reset_knots = nullptr) const override {
        return fitpack_convex_curve_c_least_squares(as<fitpack_convex_curve_c>(), smoothing, reset_knots);
    }

    /**
     * @brief set_convexity
     */
    virtual int32_t set_convexity(std::vector<double>& v) const {
        int32_t n = static_cast<int32_t>(v.size());
        return fitpack_convex_curve_c_set_convexity(as<fitpack_convex_curve_c>(), n, v.data());
    }


    int32_t comm_size() const override {
        return fitpack_convex_curve_c_comm_size(as<fitpack_convex_curve_c>());
    }

    /**
     * @brief comm_pack
     */
    void comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_convex_curve_c_comm_pack(as<fitpack_convex_curve_c>(), n, buffer.data());
    }


    /**
     * @brief comm_expand
     */
    void comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_convex_curve_c_comm_expand(as<fitpack_convex_curve_c>(), n, buffer.data());
    }


    double mse() const override {
        return fitpack_convex_curve_c_mse(as<fitpack_convex_curve_c>());
    }

    int32_t core_comm_size() const override {
        return fitpack_convex_curve_c_core_comm_size(as<fitpack_convex_curve_c>());
    }

    /**
     * @brief core_comm_pack
     */
    void core_comm_pack(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_convex_curve_c_core_comm_pack(as<fitpack_convex_curve_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    void core_comm_expand(std::vector<double>& buffer) override {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_convex_curve_c_core_comm_expand(as<fitpack_convex_curve_c>(), n, buffer.data());
    }


    void destroy_base() override {
        fitpack_convex_curve_c_destroy_base(as<fitpack_convex_curve_c>());
    }

    // ===========================================================================================
    // Component Array Accessors
    // ===========================================================================================

    /**
     * @brief Deep copy of component 'v' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> v_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_convex_curve_c_getcomp_v(as<fitpack_convex_curve_c>(), &raw, extents);
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
        fitpack_convex_curve_c_getcomp_v(as<fitpack_convex_curve_c>(), &raw, extents);
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
     * @brief Deep copy of component 'sx' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> sx_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_convex_curve_c_getcomp_sx(as<fitpack_convex_curve_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'sx'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through sx()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> sx() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_convex_curve_c_getcomp_sx(as<fitpack_convex_curve_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "sx", FX_LEN_NAME - 1);
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

    int32_t& maxtr() {
        return *fitpack_convex_curve_c_ref_maxtr(as<fitpack_convex_curve_c>());
    }
    const int32_t& maxtr() const {
        return *fitpack_convex_curve_c_ref_maxtr(as<fitpack_convex_curve_c>());
    }

    int32_t& maxbin() {
        return *fitpack_convex_curve_c_ref_maxbin(as<fitpack_convex_curve_c>());
    }
    const int32_t& maxbin() const {
        return *fitpack_convex_curve_c_ref_maxbin(as<fitpack_convex_curve_c>());
    }

protected:
    explicit fpConvexCurve(NoAlloc tag) : fpCurve(tag) {}

};

#endif /* FPCONVEXCURVE_HPP_INCLUDED */
