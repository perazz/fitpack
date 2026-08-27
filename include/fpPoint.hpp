#ifndef FPPOINT_HPP_INCLUDED
#define FPPOINT_HPP_INCLUDED

/***************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   fpPoint<dims>
!> @brief Fixed-dimension point, byte-compatible with Fortran `real(FP_REAL) :: x(dims)`
!
!   Hand-maintained (NOT emitted by the binding generator).
!
!   A `fpPoint<dims>` is a standard-layout aggregate wrapping `std::array<FP_REAL,dims>`, so a
!   `std::vector<fpPoint<dims>>` has exactly the storage of a Fortran `x(dims,m)` array in
!   column-major order: point i occupies elements [i*dims, (i+1)*dims). `fpPointData()` hands
!   that buffer straight to the rank-2 C entry points with no copy and no scatter loop.
!
!   This supersedes the pre-0.3.0 `typedef std::vector<FP_REAL> fpPoint`, which cost one heap
!   allocation per point and leaked C++ into the nominally-C `fitpack_core_c.h`.
!
!   Author: (C) Federico Perini
!> @since     26/08/2026
!
! **************************************************************************************************/

#include "fitpack_core_c.h"

#include <array>
#include <cstddef>
#include <type_traits>
#include <vector>

/**
 * @brief A point in `dims`-dimensional space.
 *
 * Aggregate-initializable with brace elision: `fpPoint<3> p{1.0, 2.0, 3.0};`
 */
template <FP_SIZE dims>
struct fpPoint
{
    static_assert(dims > 0, "fpPoint requires at least one dimension");
    static_assert(sizeof(std::array<FP_REAL, dims>) == static_cast<std::size_t>(dims) * sizeof(FP_REAL),
                  "std::array<FP_REAL,dims> must be exactly dims contiguous reals");

    //! Coordinates, in the same order as the Fortran leading dimension.
    std::array<FP_REAL, dims> coord;

    // --- Element access ---
    constexpr FP_REAL&       operator[](FP_SIZE i)       { return coord[static_cast<std::size_t>(i)]; }
    constexpr const FP_REAL& operator[](FP_SIZE i) const { return coord[static_cast<std::size_t>(i)]; }

    constexpr FP_SIZE        size() const { return dims; }
    constexpr FP_REAL*       data()       { return coord.data(); }
    constexpr const FP_REAL* data() const { return coord.data(); }

    // --- Iteration ---
    constexpr FP_REAL*       begin()        { return coord.data(); }
    constexpr FP_REAL*       end()          { return coord.data() + dims; }
    constexpr const FP_REAL* begin()  const { return coord.data(); }
    constexpr const FP_REAL* end()    const { return coord.data() + dims; }
    constexpr const FP_REAL* cbegin() const { return coord.data(); }
    constexpr const FP_REAL* cend()   const { return coord.data() + dims; }

    // --- Named accessors. Each is only instantiated when used, so z() on a 2-D point is a
    //     compile error and never a silent out-of-range read. ---
    constexpr FP_REAL& x() { static_assert(dims >= 1, "fpPoint::x() requires dims >= 1"); return coord[0]; }
    constexpr FP_REAL& y() { static_assert(dims >= 2, "fpPoint::y() requires dims >= 2"); return coord[1]; }
    constexpr FP_REAL& z() { static_assert(dims >= 3, "fpPoint::z() requires dims >= 3"); return coord[2]; }

    constexpr const FP_REAL& x() const { static_assert(dims >= 1, "fpPoint::x() requires dims >= 1"); return coord[0]; }
    constexpr const FP_REAL& y() const { static_assert(dims >= 2, "fpPoint::y() requires dims >= 2"); return coord[1]; }
    constexpr const FP_REAL& z() const { static_assert(dims >= 3, "fpPoint::z() requires dims >= 3"); return coord[2]; }
};

template <FP_SIZE dims>
constexpr bool operator==(const fpPoint<dims>& a, const fpPoint<dims>& b) { return a.coord == b.coord; }

template <FP_SIZE dims>
constexpr bool operator!=(const fpPoint<dims>& a, const fpPoint<dims>& b) { return !(a == b); }

namespace fp_detail
{
    //! Instantiated at every zero-copy call site: the layout contract that makes the cast legal.
    template <FP_SIZE dims>
    struct fpPointLayout
    {
        static_assert(sizeof(fpPoint<dims>) == static_cast<std::size_t>(dims) * sizeof(FP_REAL),
                      "fpPoint<dims> must be exactly dims contiguous FP_REALs (no padding)");
        static_assert(std::is_standard_layout<fpPoint<dims>>::value,
                      "fpPoint<dims> must be standard-layout to alias Fortran storage");
        static_assert(std::is_trivially_copyable<fpPoint<dims>>::value,
                      "fpPoint<dims> must be trivially copyable to alias Fortran storage");
        static constexpr bool value = true;
    };
}

/**
 * @brief Borrow a vector of points as a Fortran `x(dims,m)` buffer. No copy, no allocation.
 * @return `nullptr` for an empty vector, as the C entry points expect.
 */
template <FP_SIZE dims>
inline const FP_REAL* fpPointData(const std::vector<fpPoint<dims>>& points)
{
    static_assert(fp_detail::fpPointLayout<dims>::value, "fpPoint layout contract");
    return points.empty() ? nullptr : reinterpret_cast<const FP_REAL*>(points.data());
}

/**
 * @brief Flatten a vector of points into Fortran `x(dims,m)` order. One contiguous copy.
 */
template <FP_SIZE dims>
inline std::vector<FP_REAL> fpPointFlatten(const std::vector<fpPoint<dims>>& points)
{
    const FP_REAL* first = fpPointData<dims>(points);
    return std::vector<FP_REAL>(first, first + points.size() * static_cast<std::size_t>(dims));
}

/**
 * @brief Rebuild a vector of points from a Fortran `x(dims,m)` buffer.
 */
template <FP_SIZE dims>
inline std::vector<fpPoint<dims>> fpPointGather(const std::vector<FP_REAL>& flat)
{
    static_assert(fp_detail::fpPointLayout<dims>::value, "fpPoint layout contract");
    std::vector<fpPoint<dims>> points(flat.size() / static_cast<std::size_t>(dims));
    if (!points.empty())
    {
        FP_REAL* raw = reinterpret_cast<FP_REAL*>(points.data());
        for (std::size_t i = 0; i < points.size() * static_cast<std::size_t>(dims); ++i) raw[i] = flat[i];
    }
    return points;
}

// Compile-time smoke test of the layout contract for the dimensions users actually reach for.
static_assert(fp_detail::fpPointLayout<1>::value, "fpPoint<1> layout");
static_assert(fp_detail::fpPointLayout<2>::value, "fpPoint<2> layout");
static_assert(fp_detail::fpPointLayout<3>::value, "fpPoint<3> layout");
static_assert(fp_detail::fpPointLayout<4>::value, "fpPoint<4> layout");

#endif // FPPOINT_HPP_INCLUDED
