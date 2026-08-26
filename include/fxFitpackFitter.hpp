/*   ***********************************************************************************************
 *   **                                         FITPACK                                          **
 *   **                     Modern Fortran Fitting Package — C/C++ Bindings                      **
 *   ***********************************************************************************************
 *   **    fxFitpackFitter.hpp                                                                     **
 *   ** @brief Standalone C++ wrapper for fitpack_fitter (no fortran-arrays dependency)
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FXFITPACKFITTER_HPP_INCLUDED
#define FXFITPACKFITTER_HPP_INCLUDED

#include "fitpack_config.h"

#if HAVE_FXARRAY
#include "fxArrays.hpp"
#endif

#include "fitpack_fitter_c.h"
#include <string>
#include <vector>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <variant>
#include <optional>
#include <cstring>

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_fitter
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_fitter, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 */
class fxFitpackFitter {
public:
    // ===========================================================================================
    // Abstract type - cannot be instantiated directly
    // ===========================================================================================

    virtual ~fxFitpackFitter() = default;

    // ===========================================================================================
    // Utility Methods
    // ===========================================================================================

    /**
     * @brief Check if object is allocated
     */
    virtual bool is_allocated() const {
        return as<fitpack_fitter_c>()->cptr != nullptr;
    }

    /**
     * @brief Check if this is a non-owning pointer
     */
    virtual bool is_pointer() const {
        return as<fitpack_fitter_c>()->is_pointer;
    }

    /**
     * @brief Return the dynamic Fortran type name of the underlying object.
     */
    const char* fortran_type_name() const {
        const char* tag = as<fitpack_fitter_c>()->name_cptr;
        return tag ? tag : fitpack_fitter_c_typename;
    }

    /**
     * @brief Return the C wrapper struct name for this class (e.g. "fitpack_fitter_c").
     */
    virtual const char* c_type_name() const {
        return fitpack_fitter_c_c_type_name(*as<fitpack_fitter_c>());
    }

    /**
     * @brief Return the fully-qualified C++ class name for this class.
     * Owned by the C++ layer — does not round-trip through Fortran.
     */
    virtual const char* cpp_type_name() const {
        return "fxFitpackFitter";
    }

    /**
     * @brief Tag type for non-owning view constructors of polymorphic
     * concrete subtypes. Declared on roots and abstract bases so
     * descendants inherit it; variant accessors qualify the call site.
     */
    struct ViewTag {};

    /**
     * @brief Get underlying C wrapper (for interop)
     */
    fitpack_fitter_c& c_handle() {
        return *as<fitpack_fitter_c>();
    }

    /**
     * @brief Get underlying C wrapper (const, for interop)
     */
    const fitpack_fitter_c& c_handle() const {
        return *as<fitpack_fitter_c>();
    }

    // ===========================================================================================
    // Method Wrappers (standalone — no fxArray dependency)
    // ===========================================================================================

    virtual double mse() const {
        return fitpack_fitter_c_mse(as<fitpack_fitter_c>());
    }

    virtual int32_t core_comm_size() const {
        return fitpack_fitter_c_core_comm_size(as<fitpack_fitter_c>());
    }

    /**
     * @brief core_comm_pack
     */
    virtual void core_comm_pack(std::vector<double>& buffer) {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_fitter_c_core_comm_pack(as<fitpack_fitter_c>(), n, buffer.data());
    }


    /**
     * @brief core_comm_expand
     */
    virtual void core_comm_expand(std::vector<double>& buffer) {
        int32_t n = static_cast<int32_t>(buffer.size());
        fitpack_fitter_c_core_comm_expand(as<fitpack_fitter_c>(), n, buffer.data());
    }


    virtual void destroy_base() {
        fitpack_fitter_c_destroy_base(as<fitpack_fitter_c>());
    }

    virtual int32_t comm_size() const = 0;

    /**
     * @brief comm_pack (deferred)
     */
    virtual void comm_pack(std::vector<double>& buffer) = 0;


    /**
     * @brief comm_expand (deferred)
     */
    virtual void comm_expand(std::vector<double>& buffer) = 0;


    // ===========================================================================================
    // Component Array Accessors
    // ===========================================================================================

    /**
     * @brief Deep copy of component 'c' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> c_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_fitter_c_getcomp_c(as<fitpack_fitter_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'c'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through c()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> c() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_fitter_c_getcomp_c(as<fitpack_fitter_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "c", FX_LEN_NAME - 1);
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
     * @brief Deep copy of component 'wrk' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<double> wrk_vector() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_fitter_c_getcomp_wrk(as<fitpack_fitter_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const double* first = reinterpret_cast<const double*>(raw);
        return std::vector<double>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'wrk'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through wrk()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> wrk() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_fitter_c_getcomp_wrk(as<fitpack_fitter_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "wrk", FX_LEN_NAME - 1);
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
     * @brief Deep copy of component 'iwrk' as a std::vector.
     * @return An owning copy, empty when the component is unallocated.
     */
    std::vector<int32_t> iwrk_vector() const {
        int32_t* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_fitter_c_getcomp_iwrk(as<fitpack_fitter_c>(), &raw, extents);
        const int64_t total = extents[0];
        if (raw == nullptr || total <= 0) return {};
        const int32_t* first = reinterpret_cast<const int32_t*>(raw);
        return std::vector<int32_t>(first, first + total);
    }

#if HAVE_FXARRAY
    /**
     * @brief Zero-copy fxArray view of component 'iwrk'.
     *
     * The descriptor is built here from the borrowed pointer and extents, so
     * the view aliases the Fortran storage directly — nothing is copied and
     * writes through iwrk()(i) land in the object. Requires linking
     * fortran-arrays (implied by HAVE_FXARRAY, which is set from
     * __has_include("fxArrays.hpp")).
     *
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<int32_t> iwrk() const {
        int32_t* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_fitter_c_getcomp_iwrk(as<fitpack_fitter_c>(), &raw, extents);
        array_c descr = array_c_null;
        std::strncpy(descr.name, "iwrk", FX_LEN_NAME - 1);
        descr.base_address = static_cast<void*>(raw);
        descr.type = getCFITypeFlag<int32_t>();
        descr.elem_bytes = static_cast<FX_SIZE>(sizeof(int32_t));
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
        return fxArray<int32_t>(descr);
    }
#endif // HAVE_FXARRAY

    // ===========================================================================================
    // Scalar Property Accessors
    // ===========================================================================================

    int32_t& iopt() {
        return *fitpack_fitter_c_ref_iopt(as<fitpack_fitter_c>());
    }
    const int32_t& iopt() const {
        return *fitpack_fitter_c_ref_iopt(as<fitpack_fitter_c>());
    }

    double& smoothing() {
        return *fitpack_fitter_c_ref_smoothing(as<fitpack_fitter_c>());
    }
    const double& smoothing() const {
        return *fitpack_fitter_c_ref_smoothing(as<fitpack_fitter_c>());
    }

    double& fp() {
        return *fitpack_fitter_c_ref_fp(as<fitpack_fitter_c>());
    }
    const double& fp() const {
        return *fitpack_fitter_c_ref_fp(as<fitpack_fitter_c>());
    }

    int32_t& lwrk() {
        return *fitpack_fitter_c_ref_lwrk(as<fitpack_fitter_c>());
    }
    const int32_t& lwrk() const {
        return *fitpack_fitter_c_ref_lwrk(as<fitpack_fitter_c>());
    }

    int32_t& liwrk() {
        return *fitpack_fitter_c_ref_liwrk(as<fitpack_fitter_c>());
    }
    const int32_t& liwrk() const {
        return *fitpack_fitter_c_ref_liwrk(as<fitpack_fitter_c>());
    }

protected:
    /**
     * @brief Tag type for derived class constructors (no allocation)
     */
    struct NoAlloc {};

    /**
     * @brief Protected constructor for derived classes (does not allocate)
     */
    explicit fxFitpackFitter(NoAlloc) {
        // Do not allocate - derived class will allocate its own type
    }

    fitpack_fitter_c cobj = fitpack_fitter_c_null;

    template<typename T> T* as() { return static_cast<T*>(static_cast<void*>(&cobj)); }
    template<typename T> const T* as() const { return static_cast<const T*>(static_cast<const void*>(&cobj)); }

};

#endif /* FXFITPACKFITTER_HPP_INCLUDED */
