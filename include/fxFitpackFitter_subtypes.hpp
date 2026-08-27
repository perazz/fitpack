/*   ***********************************************************************************************
 *   **                                        fxArray                                          **
 *   **                                  Fortran Arrays for C++                                 **
 *   ***********************************************************************************************
 *   **    fxFitpackFitter_subtypes.hpp                                                    **
 *   ** @brief Polymorphic dispatch for class(fitpack_fitter).
 *   **
 *   ** Bundles three things, all derived from the bound concrete subtypes
 *   ** of the abstract base fxFitpackFitter:
 *   **   1. `#include`s of the abstract base + every concrete subtype
 *   **   2. `fxFitpackFitterVariant` — std::variant<Concretes...>
 *   **   3. `wrap_fitpack_fitter(handle, class_name)` — resolves the
 *   **      (base C handle, class-name string) pair into the matching
 *   **      concrete subtype variant alternative.
 *   **
 *   ** Include this header where you need to construct a polymorphic
 *   ** view — typically in parent classes that hold a class(fitpack_fitter)
 *   ** component, or in any consumer code dispatching on the dynamic type.
 *   ** The base header `fxFitpackFitter.hpp` does NOT include the
 *   ** concrete subtype headers (avoids cycles); this umbrella collects them.
 *   ***********************************************************************************************
 *   ** @author Binding Generator
 *   *********************************************************************************************** */

#ifndef FXFITPACKFITTER_SUBTYPES_HPP_INCLUDED
#define FXFITPACKFITTER_SUBTYPES_HPP_INCLUDED

#include "fxFitpackFitter.hpp"
#include "fpPeriodicCurve.hpp"
#include "fpConvexCurve.hpp"
#include "fpCurve.hpp"
#include "fpClosedCurve.hpp"
#include "fpConstrainedCurve.hpp"
#include "fpParametricCurve.hpp"
#include "fpSurface.hpp"
#include "fpGridSurface.hpp"
#include "fpGridSpline.hpp"
#include "fpParametricSurface.hpp"
#include "fpPolar.hpp"
#include "fpGridPolar.hpp"
#include "fpSphere.hpp"
#include "fpGridSphere.hpp"
#include <variant>
#include <optional>
#include <stdexcept>
#include <string>
#include <cstring>


/**
 * @brief Variant of bound concrete subtypes of class(fitpack_fitter).
 *
 * Stack-allocated tagged union — same size as one wrapper plus a tag
 * (all `*_c` wrappers share the trivial cptr+is_pointer layout). The
 * alternative-set is the C++-bound concrete subtypes enumerated at
 * generation time.
 */
using fxFitpackFitterVariant = std::variant<
    fpPeriodicCurve,
    fpConvexCurve,
    fpCurve,
    fpClosedCurve,
    fpConstrainedCurve,
    fpParametricCurve,
    fpSurface,
    fpGridSurface,
    fpGridSpline,
    fpParametricSurface,
    fpPolar,
    fpGridPolar,
    fpSphere,
    fpGridSphere>;

/**
 * @brief Resolve a (handle, class_name) pair into the right concrete
 * variant alternative.
 *
 * @param handle      A base-class C handle (e.g. from a polymorphic
 *                    component accessor — any caller that has a
 *                    `fitpack_fitter_c` plus the dynamic type name).
 * @param class_name  The dynamic Fortran type name string returned
 *                    alongside the handle (e.g. "fitpack_periodic_curve").
 *
 * @return std::nullopt when @p handle.cptr is null (component
 *         unallocated, or `class default` arm fired in Fortran).
 *
 * @throws std::runtime_error if @p class_name is null but the handle
 *         isn't (a generator bug — the polymorphic accessor's
 *         select_type fanout missed an arm) or if @p class_name does
 *         not match any registered subtype (binding/YAML drift —
 *         regenerate bindings).
 *
 * Currently constructs a non-owning view of the handle's storage
 * (is_pointer=true on the inner C handle). A future overload may
 * construct an owned variant from a transferred handle.
 */
inline std::optional<fxFitpackFitterVariant>
wrap_fitpack_fitter(fitpack_fitter_c handle, const char* class_name) {
    if (handle.cptr == nullptr) return std::nullopt;
    if (class_name == nullptr) {
        throw std::runtime_error(
            "wrap_fitpack_fitter: handle non-null but class_name null "
            "(generator bug: select_type fanout missed an arm).");
    }
    if (std::strcmp(class_name, "fitpack_periodic_curve") == 0)
        return fxFitpackFitterVariant(fpPeriodicCurve(
            *reinterpret_cast<fitpack_periodic_curve_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_convex_curve") == 0)
        return fxFitpackFitterVariant(fpConvexCurve(
            *reinterpret_cast<fitpack_convex_curve_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_curve") == 0)
        return fxFitpackFitterVariant(fpCurve(
            *reinterpret_cast<fitpack_curve_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_closed_curve") == 0)
        return fxFitpackFitterVariant(fpClosedCurve(
            *reinterpret_cast<fitpack_closed_curve_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_constrained_curve") == 0)
        return fxFitpackFitterVariant(fpConstrainedCurve(
            *reinterpret_cast<fitpack_constrained_curve_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_parametric_curve") == 0)
        return fxFitpackFitterVariant(fpParametricCurve(
            *reinterpret_cast<fitpack_parametric_curve_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_surface") == 0)
        return fxFitpackFitterVariant(fpSurface(
            *reinterpret_cast<fitpack_surface_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_grid_surface") == 0)
        return fxFitpackFitterVariant(fpGridSurface(
            *reinterpret_cast<fitpack_grid_surface_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_gridded_spline") == 0)
        return fxFitpackFitterVariant(fpGridSpline(
            *reinterpret_cast<fitpack_gridded_spline_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_parametric_surface") == 0)
        return fxFitpackFitterVariant(fpParametricSurface(
            *reinterpret_cast<fitpack_parametric_surface_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_polar") == 0)
        return fxFitpackFitterVariant(fpPolar(
            *reinterpret_cast<fitpack_polar_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_grid_polar") == 0)
        return fxFitpackFitterVariant(fpGridPolar(
            *reinterpret_cast<fitpack_grid_polar_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_sphere") == 0)
        return fxFitpackFitterVariant(fpSphere(
            *reinterpret_cast<fitpack_sphere_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    if (std::strcmp(class_name, "fitpack_grid_sphere") == 0)
        return fxFitpackFitterVariant(fpGridSphere(
            *reinterpret_cast<fitpack_grid_sphere_c*>(&handle),
            fxFitpackFitter::ViewTag{}));
    throw std::runtime_error(
        std::string("wrap_fitpack_fitter: unbound dynamic type '") +
        class_name + "' (regenerate bindings).");
}


#endif /* FXFITPACKFITTER_SUBTYPES_HPP_INCLUDED */
