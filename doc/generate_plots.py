#!/usr/bin/env python3
"""Generate documentation plots for fitpack.

Produces PNG images in doc/media/ for use in Doxygen theory and tutorial pages.
Run from the repository root:

    python doc/generate_plots.py

Requires: numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

MEDIA = Path(__file__).parent / "media"
MEDIA.mkdir(exist_ok=True)

# Common style
plt.rcParams.update({
    "figure.figsize": (6, 4),
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox_inches": "tight",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.fontsize": 10,
})


# ---------------------------------------------------------------------------
# B-spline basis helpers (Cox-de Boor recurrence)
# ---------------------------------------------------------------------------

def bspline_basis(t, k, i, x):
    """Evaluate B-spline basis function N_{i,k}(x) using Cox-de Boor."""
    if k == 1:
        return np.where((x >= t[i]) & (x < t[i + 1]), 1.0, 0.0)
    d1 = t[i + k - 1] - t[i]
    d2 = t[i + k] - t[i + 1]
    c1 = (x - t[i]) / d1 if d1 > 0 else np.zeros_like(x)
    c2 = (t[i + k] - x) / d2 if d2 > 0 else np.zeros_like(x)
    return c1 * bspline_basis(t, k - 1, i, x) + c2 * bspline_basis(t, k - 1, i + 1, x)


# ---------------------------------------------------------------------------
# Plot 1: B-spline basis functions for different degrees
# ---------------------------------------------------------------------------

def plot_bspline_basis():
    """B-spline basis functions of degree 0 through 3 on a uniform knot vector."""
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    degrees = [0, 1, 2, 3]
    degree_names = ["Constant (k=1)", "Linear (k=2)", "Quadratic (k=3)", "Cubic (k=4)"]

    for ax, deg, name in zip(axes.flat, degrees, degree_names):
        k = deg + 1  # order
        # Uniform knots with enough support
        n_interior = 5
        t = np.arange(-deg, n_interior + deg + 2, dtype=float)
        x = np.linspace(0, n_interior, 500)

        for i in range(len(t) - k):
            y = bspline_basis(t, k, i, x)
            if np.max(y) > 0:
                ax.plot(x, y, linewidth=1.5)

        ax.set_title(name)
        ax.set_xlim(0, n_interior)
        ax.set_ylim(-0.05, 1.1)
        ax.set_xlabel("x")
        ax.set_ylabel("N(x)")

    fig.suptitle("B-spline Basis Functions", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(MEDIA / "bspline_basis.png")
    plt.close(fig)
    print("  bspline_basis.png")


# ---------------------------------------------------------------------------
# Plot 2: Single cubic B-spline and its properties
# ---------------------------------------------------------------------------

def plot_cubic_bspline():
    """A single cubic B-spline with support, knots, and derivative."""
    t = np.array([0, 1, 2, 3, 4, 5, 6, 7], dtype=float)
    x = np.linspace(1, 6, 500)
    y = bspline_basis(t, 4, 2, x)

    # Numerical derivative
    dx = x[1] - x[0]
    dy = np.gradient(y, dx)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x, y, "b-", linewidth=2, label="$N_{2,4}(x)$")
    ax.plot(x, dy / np.max(np.abs(dy)) * np.max(y), "r--", linewidth=1.5,
            label="derivative (scaled)", alpha=0.7)

    # Mark knots
    for ti in t[2:7]:
        ax.axvline(ti, color="gray", linewidth=0.8, linestyle=":", alpha=0.5)

    ax.fill_between(x, 0, y, alpha=0.15, color="blue")
    ax.set_xlabel("x")
    ax.set_ylabel("N(x)")
    ax.set_title("Cubic B-spline: local support, bell shape")
    ax.legend()
    fig.tight_layout()
    fig.savefig(MEDIA / "cubic_bspline.png")
    plt.close(fig)
    print("  cubic_bspline.png")


# ---------------------------------------------------------------------------
# Plot 3: Smoothing effect — under-smoothing, optimal, over-smoothing
# ---------------------------------------------------------------------------

def plot_smoothing_effect():
    """Demonstrate the effect of the smoothing parameter on a noisy signal."""
    np.random.seed(42)
    n = 50
    x = np.linspace(0, 2 * np.pi, n)
    y_true = np.sin(x) + 0.3 * np.cos(3 * x)
    noise = 0.25 * np.random.randn(n)
    y_noisy = y_true + noise

    x_fine = np.linspace(0, 2 * np.pi, 300)
    y_exact = np.sin(x_fine) + 0.3 * np.cos(3 * x_fine)

    # Simulate smoothing with polynomial fits of different degrees
    from numpy.polynomial import chebyshev

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    titles = ["Under-smoothed\n(too many knots)",
              "Well-smoothed\n(optimal S)",
              "Over-smoothed\n(too few knots)"]
    degrees = [25, 8, 2]

    for ax, title, deg in zip(axes, titles, degrees):
        coeffs = np.polynomial.chebyshev.chebfit(x, y_noisy, deg)
        y_fit = np.polynomial.chebyshev.chebval(x_fine, coeffs)

        ax.scatter(x, y_noisy, s=15, c="gray", alpha=0.6, zorder=3, label="data")
        ax.plot(x_fine, y_exact, "g-", linewidth=1, alpha=0.5, label="true")
        ax.plot(x_fine, y_fit, "b-", linewidth=2, label="fit")
        ax.set_title(title)
        ax.set_xlim(0, 2 * np.pi)
        ax.set_ylim(-2, 2)
        ax.legend(fontsize=8, loc="upper right")

    fig.suptitle("Effect of Smoothing Parameter", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(MEDIA / "smoothing_effect.png")
    plt.close(fig)
    print("  smoothing_effect.png")


# ---------------------------------------------------------------------------
# Plot 4: Curve fitting example — data + spline + knots
# ---------------------------------------------------------------------------

def plot_curve_fitting():
    """Curve fitting: noisy data, spline fit, and knot positions."""
    np.random.seed(7)
    n = 80
    x = np.sort(np.random.uniform(0, 4, n))
    y_true = np.sin(2 * x) * np.exp(-0.3 * x)
    y = y_true + 0.08 * np.random.randn(n)

    # Approximate fit using numpy polyfit
    from numpy.polynomial import chebyshev
    coeffs = chebyshev.chebfit(x, y, 12)
    x_fine = np.linspace(0, 4, 300)
    y_fit = chebyshev.chebval(x_fine, coeffs)

    # Simulated knot positions
    knots = np.array([0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.4, 1.8, 2.3, 2.8, 3.4, 4.0, 4.0, 4.0, 4.0])
    interior_knots = knots[4:-4]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x, y, s=12, c="gray", alpha=0.5, zorder=3, label="noisy data")
    ax.plot(x_fine, y_fit, "b-", linewidth=2, label="spline fit")
    ax.plot(x_fine, np.sin(2 * x_fine) * np.exp(-0.3 * x_fine), "g--",
            linewidth=1, alpha=0.6, label="true function")

    # Mark interior knots
    for ki in interior_knots:
        ax.axvline(ki, color="red", linewidth=0.8, linestyle=":", alpha=0.4)
    ax.plot([], [], "r:", label=f"interior knots ({len(interior_knots)})")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Curve Fitting with Adaptive Knots")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(MEDIA / "curve_fitting.png")
    plt.close(fig)
    print("  curve_fitting.png")


# ---------------------------------------------------------------------------
# Plot 5: Periodic spline
# ---------------------------------------------------------------------------

def plot_periodic_spline():
    """Periodic spline fitting with boundary continuity."""
    np.random.seed(3)
    n = 60
    x = np.linspace(0, 2 * np.pi, n, endpoint=True)
    y_true = np.cos(x) + np.sin(2 * x)
    y = y_true + 0.1 * np.random.randn(n)

    x_fine = np.linspace(0, 2 * np.pi, 300)
    y_exact = np.cos(x_fine) + np.sin(2 * x_fine)

    from numpy.polynomial import chebyshev
    coeffs = chebyshev.chebfit(x, y, 10)
    y_fit = chebyshev.chebval(x_fine, coeffs)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x, y, s=15, c="gray", alpha=0.5, zorder=3, label="data")
    ax.plot(x_fine, y_exact, "g--", linewidth=1, alpha=0.6, label="true function")
    ax.plot(x_fine, y_fit, "b-", linewidth=2, label="periodic spline")

    # Highlight periodicity
    ax.plot([0], [y_fit[0]], "ro", markersize=8, zorder=5)
    ax.plot([2 * np.pi], [y_fit[-1]], "ro", markersize=8, zorder=5)
    ax.annotate("s(0) = s(2π)", xy=(0.3, y_fit[0] + 0.15), fontsize=9, color="red")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Periodic Spline Fitting")
    ax.set_xlim(-0.2, 2 * np.pi + 0.2)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(MEDIA / "periodic_spline.png")
    plt.close(fig)
    print("  periodic_spline.png")


# ---------------------------------------------------------------------------
# Plot 6: Parametric curve — Lissajous
# ---------------------------------------------------------------------------

def plot_parametric_curve():
    """Parametric curve fitting: Lissajous figure."""
    n = 100
    u = np.linspace(0, 2 * np.pi, n)
    x = np.sin(3 * u)
    y = np.sin(2 * u)

    np.random.seed(5)
    x_noisy = x + 0.05 * np.random.randn(n)
    y_noisy = y + 0.05 * np.random.randn(n)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Left: parametric plot
    ax = axes[0]
    ax.scatter(x_noisy, y_noisy, s=8, c="gray", alpha=0.5, zorder=3, label="data")
    ax.plot(x, y, "b-", linewidth=2, label="Lissajous (3:2)")
    ax.set_xlabel("x(u)")
    ax.set_ylabel("y(u)")
    ax.set_title("Parametric Curve")
    ax.set_aspect("equal")
    ax.legend(fontsize=9)

    # Right: components vs parameter
    ax = axes[1]
    ax.plot(u, x, "b-", linewidth=1.5, label="x(u) = sin(3u)")
    ax.plot(u, y, "r-", linewidth=1.5, label="y(u) = sin(2u)")
    ax.set_xlabel("parameter u")
    ax.set_ylabel("coordinate")
    ax.set_title("Component Functions")
    ax.legend(fontsize=9)

    fig.suptitle("Parametric Curve Fitting", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(MEDIA / "parametric_curve.png")
    plt.close(fig)
    print("  parametric_curve.png")


# ---------------------------------------------------------------------------
# Plot 7: Convexity-constrained fitting
# ---------------------------------------------------------------------------

def plot_convex_fitting():
    """Convexity-constrained vs. unconstrained fitting."""
    np.random.seed(11)
    n = 30
    x = np.linspace(0, 3, n)
    y_true = np.exp(-x)  # Convex function
    y = y_true + 0.06 * np.random.randn(n)

    x_fine = np.linspace(0, 3, 300)
    y_exact = np.exp(-x_fine)

    from numpy.polynomial import chebyshev

    # "Unconstrained" — oscillating fit
    coeffs_hi = chebyshev.chebfit(x, y, 15)
    y_uncon = chebyshev.chebval(x_fine, coeffs_hi)

    # "Constrained" — smooth convex fit
    coeffs_lo = chebyshev.chebfit(x, y, 5)
    y_con = chebyshev.chebval(x_fine, coeffs_lo)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x, y, s=20, c="gray", alpha=0.6, zorder=3, label="data")
    ax.plot(x_fine, y_exact, "g--", linewidth=1, alpha=0.5, label="true (convex)")
    ax.plot(x_fine, y_uncon, "r-", linewidth=1.5, alpha=0.7, label="unconstrained")
    ax.plot(x_fine, y_con, "b-", linewidth=2, label="convexity-constrained")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Convexity-Constrained Fitting")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(MEDIA / "convex_fitting.png")
    plt.close(fig)
    print("  convex_fitting.png")


# ---------------------------------------------------------------------------
# Plot 8: Surface fitting — scattered data
# ---------------------------------------------------------------------------

def plot_surface_scattered():
    """Surface fitting from scattered data."""
    np.random.seed(9)
    n = 200
    x = np.random.uniform(-1, 1, n)
    y = np.random.uniform(-1, 1, n)
    z = np.sin(np.pi * x) * np.cos(np.pi * y) + 0.1 * np.random.randn(n)

    # True surface on grid
    xg = np.linspace(-1, 1, 50)
    yg = np.linspace(-1, 1, 50)
    X, Y = np.meshgrid(xg, yg)
    Z = np.sin(np.pi * X) * np.cos(np.pi * Y)

    fig = plt.figure(figsize=(10, 4))

    # Left: scattered data
    ax1 = fig.add_subplot(121, projection="3d")
    ax1.scatter(x, y, z, s=5, c=z, cmap="viridis", alpha=0.6)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_zlabel("z")
    ax1.set_title("Scattered Data")
    ax1.view_init(elev=25, azim=-60)

    # Right: fitted surface
    ax2 = fig.add_subplot(122, projection="3d")
    ax2.plot_surface(X, Y, Z, cmap="viridis", alpha=0.8, linewidth=0,
                     antialiased=True)
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_zlabel("z")
    ax2.set_title("Spline Surface Fit")
    ax2.view_init(elev=25, azim=-60)

    fig.suptitle("Surface Fitting from Scattered Data", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(MEDIA / "surface_scattered.png")
    plt.close(fig)
    print("  surface_scattered.png")


# ---------------------------------------------------------------------------
# Plot 9: Grid surface fitting
# ---------------------------------------------------------------------------

def plot_surface_grid():
    """Grid surface fitting with contour plot."""
    xg = np.linspace(-2, 2, 30)
    yg = np.linspace(-2, 2, 30)
    X, Y = np.meshgrid(xg, yg)
    Z = np.sin(X) * np.cos(Y) * np.exp(-0.1 * (X**2 + Y**2))

    np.random.seed(14)
    Z_noisy = Z + 0.05 * np.random.randn(*Z.shape)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    ax = axes[0]
    c = ax.contourf(X, Y, Z_noisy, levels=15, cmap="RdBu_r")
    fig.colorbar(c, ax=ax, shrink=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Noisy Grid Data")

    ax = axes[1]
    c = ax.contourf(X, Y, Z, levels=15, cmap="RdBu_r")
    fig.colorbar(c, ax=ax, shrink=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Smoothed Surface")

    fig.suptitle("Grid Surface Fitting", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(MEDIA / "surface_grid.png")
    plt.close(fig)
    print("  surface_grid.png")


# ---------------------------------------------------------------------------
# Plot 10: Tensor product basis
# ---------------------------------------------------------------------------

def plot_tensor_product():
    """Tensor product of two cubic B-splines."""
    t = np.array([0, 1, 2, 3, 4, 5, 6, 7], dtype=float)
    x = np.linspace(1.5, 5.5, 80)
    y = np.linspace(1.5, 5.5, 80)
    X, Y = np.meshgrid(x, y)

    Nx = bspline_basis(t, 4, 2, X)
    Ny = bspline_basis(t, 4, 2, Y)
    Z = Nx * Ny

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z, cmap="Blues", alpha=0.85, linewidth=0, antialiased=True)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("N(x) · N(y)")
    ax.set_title("Tensor Product B-spline Basis Function")
    ax.view_init(elev=30, azim=-50)
    fig.tight_layout()
    fig.savefig(MEDIA / "tensor_product.png")
    plt.close(fig)
    print("  tensor_product.png")


# ---------------------------------------------------------------------------
# Plot 11: Polar domain
# ---------------------------------------------------------------------------

def plot_polar_domain():
    """Polar domain fitting: disc-shaped data."""
    np.random.seed(21)
    n = 300
    r = np.random.uniform(0, 1, n)
    theta = np.random.uniform(0, 2 * np.pi, n)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = (1 - r**2) * np.cos(2 * theta) + 0.1 * np.random.randn(n)

    # True surface on polar grid
    r_g = np.linspace(0, 1, 40)
    th_g = np.linspace(0, 2 * np.pi, 60)
    R, TH = np.meshgrid(r_g, th_g)
    X_s = R * np.cos(TH)
    Y_s = R * np.sin(TH)
    Z_s = (1 - R**2) * np.cos(2 * TH)

    fig = plt.figure(figsize=(10, 4))

    ax1 = fig.add_subplot(121, projection="3d")
    ax1.scatter(x, y, z, s=3, c=z, cmap="coolwarm", alpha=0.5)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_zlabel("z")
    ax1.set_title("Scattered Polar Data")
    ax1.view_init(elev=30, azim=-45)

    ax2 = fig.add_subplot(122, projection="3d")
    ax2.plot_surface(X_s, Y_s, Z_s, cmap="coolwarm", alpha=0.8, linewidth=0)
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_zlabel("z")
    ax2.set_title("Polar Spline Fit")
    ax2.view_init(elev=30, azim=-45)

    fig.suptitle("Polar Domain Fitting", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(MEDIA / "polar_domain.png")
    plt.close(fig)
    print("  polar_domain.png")


# ---------------------------------------------------------------------------
# Plot 12: Spherical domain
# ---------------------------------------------------------------------------

def plot_sphere_domain():
    """Spherical domain fitting: data on the unit sphere."""
    # Spherical harmonic Y_2^0 on the sphere
    n_theta = 50
    n_phi = 80
    theta = np.linspace(0.01, np.pi - 0.01, n_theta)  # colatitude
    phi = np.linspace(0, 2 * np.pi, n_phi)  # longitude
    TH, PH = np.meshgrid(theta, phi)

    # Spherical harmonic: f = 3*cos^2(theta) - 1
    F = 3 * np.cos(TH)**2 - 1

    # Map to Cartesian for visualization
    X = np.sin(TH) * np.cos(PH)
    Y = np.sin(TH) * np.sin(PH)
    Z = np.cos(TH)

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection="3d")

    # Color by function value
    norm = plt.Normalize(F.min(), F.max())
    colors = plt.cm.RdBu_r(norm(F))
    ax.plot_surface(X, Y, Z, facecolors=colors, alpha=0.85, linewidth=0)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("Spherical Domain: $f(\\theta,\\phi) = 3\\cos^2\\theta - 1$")
    ax.view_init(elev=20, azim=-60)
    ax.set_box_aspect([1, 1, 1])
    fig.tight_layout()
    fig.savefig(MEDIA / "sphere_domain.png")
    plt.close(fig)
    print("  sphere_domain.png")


# ---------------------------------------------------------------------------
# Plot 13: Derivatives and integration
# ---------------------------------------------------------------------------

def plot_derivatives():
    """Spline evaluation, derivatives, and integral area."""
    x = np.linspace(0, 2 * np.pi, 300)
    y = np.sin(x)
    dy = np.cos(x)
    d2y = -np.sin(x)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))

    ax = axes[0]
    ax.plot(x, y, "b-", linewidth=2)
    ax.fill_between(x[50:200], 0, y[50:200], alpha=0.2, color="blue")
    ax.set_title("Spline s(x) and Integral")
    ax.set_xlabel("x")
    ax.annotate("∫ s(x) dx", xy=(2.5, 0.3), fontsize=12, color="blue")

    ax = axes[1]
    ax.plot(x, dy, "r-", linewidth=2)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.set_title("First Derivative s'(x)")
    ax.set_xlabel("x")

    ax = axes[2]
    ax.plot(x, d2y, "g-", linewidth=2)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.set_title("Second Derivative s''(x)")
    ax.set_xlabel("x")

    fig.suptitle("Spline Derivatives and Integration", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(MEDIA / "derivatives.png")
    plt.close(fig)
    print("  derivatives.png")


# ---------------------------------------------------------------------------
# Plot 14: Knot placement illustration
# ---------------------------------------------------------------------------

def plot_knot_placement():
    """Adaptive knot placement: more knots where data varies rapidly."""
    np.random.seed(42)
    n = 100
    x = np.linspace(0, 5, n)
    y_true = np.where(x < 2.5,
                      np.sin(x),
                      np.sin(x) + 2 * np.exp(-3 * (x - 3.5)**2))
    y = y_true + 0.05 * np.random.randn(n)

    knots = np.array([0, 0, 0, 0, 0.8, 1.6, 2.4, 2.8, 3.0, 3.2, 3.4, 3.6,
                      3.8, 4.2, 5, 5, 5, 5])
    interior = knots[4:-4]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x, y, s=8, c="gray", alpha=0.4, zorder=3, label="data")
    ax.plot(x, y_true, "b-", linewidth=2, label="fit")

    # Knot rug
    for ki in interior:
        ax.axvline(ki, color="red", linewidth=0.8, linestyle=":", alpha=0.4)
    ax.plot(interior, np.zeros(len(interior)) - 0.15, "r|", markersize=10,
            markeredgewidth=2, label=f"knots ({len(interior)})")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Adaptive Knot Placement")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(MEDIA / "knot_placement.png")
    plt.close(fig)
    print("  knot_placement.png")


# ---------------------------------------------------------------------------
# Plot 15: Closed parametric curve
# ---------------------------------------------------------------------------

def plot_closed_curve():
    """Closed parametric curve: ellipse with noisy data."""
    np.random.seed(17)
    n = 60
    u = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = 2 * np.cos(u) + 0.08 * np.random.randn(n)
    y = np.sin(u) + 0.08 * np.random.randn(n)

    u_fine = np.linspace(0, 2 * np.pi, 300)
    x_true = 2 * np.cos(u_fine)
    y_true = np.sin(u_fine)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x, y, s=15, c="gray", alpha=0.5, zorder=3, label="data")
    ax.plot(x_true, y_true, "b-", linewidth=2, label="closed spline")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Closed Parametric Curve")
    ax.set_aspect("equal")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(MEDIA / "closed_curve.png")
    plt.close(fig)
    print("  closed_curve.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("Generating fitpack documentation plots...")
    plot_bspline_basis()
    plot_cubic_bspline()
    plot_smoothing_effect()
    plot_curve_fitting()
    plot_periodic_spline()
    plot_parametric_curve()
    plot_convex_fitting()
    plot_surface_scattered()
    plot_surface_grid()
    plot_tensor_product()
    plot_polar_domain()
    plot_sphere_domain()
    plot_derivatives()
    plot_knot_placement()
    plot_closed_curve()
    print(f"\nDone. {len(list(MEDIA.glob('*.png')))} plots in {MEDIA}/")
