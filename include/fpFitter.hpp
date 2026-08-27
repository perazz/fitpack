/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fpFitter.hpp (class fpFitter)
!> @brief Standalone C++ wrapper for fitpack_fitter (no fortran-arrays dependency)
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

#ifndef FPFITTER_HPP_INCLUDED
#define FPFITTER_HPP_INCLUDED

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

/**
 * @brief Standalone C++ RAII wrapper for Fortran fitpack_fitter
 *
 * This class provides automatic memory management (RAII) and a natural C++ API
 * for the underlying Fortran fitpack_fitter, without requiring the fortran-arrays library.
 * Array arguments use std::vector<T> for standalone operation.
 */
class fpFitter {
public:
    // ===========================================================================================
    // Abstract type - cannot be instantiated directly
    // ===========================================================================================

    virtual ~fpFitter() = default;

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
        return "fpFitter";
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
    virtual void core_comm_pack(std::vector<double>& buffer) const {
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
    virtual void comm_pack(std::vector<double>& buffer) const = 0;


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
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through c()(i) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> c() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_fitter_c_getcomp_c(as<fitpack_fitter_c>(), &raw, extents);
        FX_SIZE bounds[2];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 1; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "c", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(1), bounds);
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
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through wrk()(i) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<double> wrk() const {
        double* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_fitter_c_getcomp_wrk(as<fitpack_fitter_c>(), &raw, extents);
        FX_SIZE bounds[2];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 1; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "wrk", static_cast<void*>(raw),
                         getCFITypeFlag<double>(),
                         static_cast<FX_SIZE>(sizeof(double)),
                         static_cast<FX_RANK>(1), bounds);
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
     * `array_c_from_ptr` builds the descriptor from the borrowed pointer and
     * bounds, so the view aliases the Fortran storage directly — nothing is
     * copied and writes through iwrk()(i) land in the object.
     *
     * @note Only available when the Fortran-Arrays library is present and
     *       enabled (HAVE_FXARRAY).
     * @warning The view is invalidated by any refit, assignment or destroy
     *          call on this object; re-read it rather than retaining it.
     */
    fxArray<int32_t> iwrk() const {
        int32_t* raw = nullptr;
        int64_t extents[1] = {0};
        fitpack_fitter_c_getcomp_iwrk(as<fitpack_fitter_c>(), &raw, extents);
        FX_SIZE bounds[2];   // (lower, upper) per dimension, Fortran lbound 1
        for (int k = 0; k < 1; ++k) {
            bounds[2 * k] = 1;
            bounds[2 * k + 1] = static_cast<FX_SIZE>(extents[k]);
        }
        array_c descr = array_c_null;
        array_c_from_ptr(&descr, "iwrk", static_cast<void*>(raw),
                         getCFITypeFlag<int32_t>(),
                         static_cast<FX_SIZE>(sizeof(int32_t)),
                         static_cast<FX_RANK>(1), bounds);
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
    explicit fpFitter(NoAlloc) {
        // Do not allocate - derived class will allocate its own type
    }

    fitpack_fitter_c cobj = fitpack_fitter_c_null;

    template<typename T> T* as() { return static_cast<T*>(static_cast<void*>(&cobj)); }
    template<typename T> const T* as() const { return static_cast<const T*>(static_cast<const void*>(&cobj)); }

};

#endif /* FPFITTER_HPP_INCLUDED */
