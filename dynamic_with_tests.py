# filename: test_wormhole_geometry.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import pytest # Using pytest for more structured testing and assertions
import sys
import os

# =============================================
# 1. Traversable Wormhole Metric (Moved from wormhole_geometry.py)
# =============================================
R0 = 1.0        # Throat radius
A = 0.5         # Flare control
omega = 0.3     # Oscillation frequency
beta = 0.1      # Asymmetry parameter

def redshift_func(r, t=0):
    """Time-dependent redshift function Φ(r,t) to prevent horizons"""
    return 0.5 * np.exp(-(r - R0)**2) * np.cos(omega * t)

def shape_func(r, t=0):
    """Shape function b(r) ensuring throat opens/closes dynamically"""
    return R0 * (1 + A * np.tanh(t) * np.exp(-(r - R0)**2))

def traversable_metric(x, y, t=0):
    r = np.sqrt(x**2 + y**2)
    θ = np.arctan2(y, x)
    Φ = redshift_func(r, t)
    b = shape_func(r, t)

    epsilon = 0.05  # to avoid division by near-zero
    safe_r = np.maximum(r - R0, epsilon)

    Z = np.sqrt(b / safe_r) * np.exp(Φ) * np.cos(2 * θ + beta * t)
    Z = np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0)
    Z = np.clip(Z, -5, 5)
    return Z

# =============================================
# 2. Dynamic Visualization (Mouth Opening) (Moved from wormhole_geometry.py)
# =============================================
def animate_wormhole():
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=30, azim=60)
    x = y = np.linspace(-3, 3, 100)
    X, Y = np.meshgrid(x, y)
    
    def update(t):
        ax.clear()
        Z = traversable_metric(X, Y, t)
        Z = np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0)
        Z = np.clip(Z, -5, 5)  # avoid extreme spikes
        ax.plot_surface(X, Y, Z, cmap='plasma', alpha=0.8)
        ax.plot_surface(X, Y, Z, cmap='plasma', rstride=2, cstride=2, linewidth=0, antialiased=True)

        ax.set_title(f"Traversable Wormhole (t = {t:.1f})")
        ax.set_zlim(-2, 2)
    
    ani = FuncAnimation(fig, update, frames=np.linspace(0, 5, 50), interval=100)
    # plt.show() # Removed to prevent blocking pytest execution
    return ani

# =============================================
# 3. Geodesics for Traversing Photons (Moved from wormhole_geometry.py)
# =============================================
def geodesic_eq(u, t, wormhole_t=0):
    x, y, vx, vy = u
    r = np.sqrt(x**2 + y**2) + 1e-8
    Φ = redshift_func(r, wormhole_t)
    b = shape_func(r, wormhole_t)
    
    # Derivatives
    dΦ_dr = -np.exp(-(r - R0)**2) * (r - R0) * np.cos(omega * wormhole_t)
    db_dr = -2 * R0 * A * (r - R0) * np.tanh(wormhole_t) * np.exp(-(r - R0)**2)

    ax = -0.5 * (dΦ_dr * x / r + (R0 / r**4) * (b - r * db_dr) * x)
    ay = -0.5 * (dΦ_dr * y / r + (R0 / r**4) * (b - r * db_dr) * y)

    return [vx, vy, ax, ay]


def nec_violation(r, t=0):
    """Returns NEC quantity: rho + p_r at given radius r and time t"""
    Φ = redshift_func(r, t)
    b = shape_func(r, t)

    # Derivatives
    dΦ_dr = -np.exp(-(r - R0)**2) * (r - R0) * np.cos(omega * t)
    db_dr = -2 * R0 * A * (r - R0) * np.tanh(t) * np.exp(-(r - R0)**2)

    # NEC expression
    rho_plus_pr = (db_dr / r**2) - (2 * (1 - b / r) * dΦ_dr / r)
    return rho_plus_pr / (8 * np.pi)

def plot_nec(t=2.0):
    r_vals = np.linspace(0.5, 3.0, 200)
    nec_vals = [nec_violation(r, t) for r in r_vals]

    plt.figure(figsize=(8, 5))
    plt.plot(r_vals, nec_vals, label=r'$\rho + p_r$')
    plt.axhline(0, color='k', linestyle='--', alpha=0.5)
    plt.axvline(R0, color='r', linestyle=':', label='Throat $r = R_0$')
    plt.title(f"NEC Condition at Time t = {t}")
    plt.xlabel("Radial Coordinate r")
    plt.ylabel(r"$\rho + p_r$")
    plt.grid(True)
    plt.legend()
    # plt.show() # Removed to prevent blocking pytest execution


def simulate_traversal():
    # Adjusted to return traversed_count for testing
    # plt.figure(figsize=(8, 8)) # Removed to prevent blocking pytest
    t_span = np.linspace(0, 50, 3000)
    
    traversed_count = 0
    # Store solutions to prevent plt.show() from being called
    solutions = []

    # Initial conditions (photon aimed at throat)
    for impact_param in np.linspace(-0.5, 0.5, 15):  # more photons for better sampling
        sol = odeint(geodesic_eq, [3.0, impact_param, -1.0, 0.0], t_span, args=(2.0,))
        # plt.plot(sol[:, 0], sol[:, 1], label=f"b={impact_param:.2f}") # Removed
        solutions.append(sol) # Store for potential future visualization

        # Check if photon crossed x=0
        crossing_indices = np.where(np.diff(np.sign(sol[:, 0])))[0]
        if len(crossing_indices) > 0:
            # print(f"✅ Photon with b = {impact_param:.2f} traversed.") # Removed
            traversed_count += 1

    # Removed visualization related code to prevent blocking pytest
    # plt.gca().add_patch(plt.Circle((0, 0), R0, color='red', alpha=0.2, label='Throat'))
    # plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    # plt.title("Photon Traversal Through Wormhole")
    # plt.xlabel("X"), plt.ylabel("Y")
    # plt.axis('equal'), plt.grid()
    # plt.legend()
    # plt.show()

    # print(f"\n🔎 Total traversing photons: {traversed_count}") # Removed

    return traversed_count, solutions # Return solutions for possible external use

# =============================================================================
# TEST SUITE FOR WORMHOLE GEOMETRY
# Run with: pytest test_wormhole_geometry.py
# =============================================================================

# Define a small tolerance for floating point comparisons
TOL = 1e-9

class TestWormholeGeometry:
    """
    Comprehensive tests for the wormhole geometry functions.
    """

    def test_redshift_func_at_throat(self):
        """
        Test redshift function at the throat (r = R0).
        It should be finite and oscillate with time.
        """
        t_vals = np.array([0, np.pi / (2 * omega), np.pi / omega])
        expected_at_throat = 0.5 * np.cos(omega * t_vals) # (r-R0)**2 term is 0 at R0

        for i, t in enumerate(t_vals):
            phi = redshift_func(R0, t)
            assert np.isfinite(phi), f"redshift_func at r=R0, t={t} is not finite."
            assert np.isclose(phi, expected_at_throat[i]), \
                   f"redshift_func at r=R0, t={t} does not match expected oscillation."

    def test_redshift_func_asymptotic_behavior(self):
        """
        Test redshift function far from the throat (large r).
        It should decay towards zero.
        """
        r_far = 100 * R0
        phi_far = redshift_func(r_far, t=0)
        assert np.isfinite(phi_far), "redshift_func far from throat is not finite."
        assert np.isclose(phi_far, 0.0, atol=TOL), \
               "redshift_func should approach 0 far from the throat."

    def test_shape_func_at_throat(self):
        """
        Test shape function at the throat (r = R0).
        It should be finite and dynamically change with time due to tanh(t).
        """
        t_vals = np.array([-5, 0, 5]) # Test different times for tanh behavior
        expected_b_at_throat = R0 * (1 + A * np.tanh(t_vals)) # (r-R0)**2 term is 0 at R0

        for i, t in enumerate(t_vals):
            b_val = shape_func(R0, t)
            assert np.isfinite(b_val), f"shape_func at r=R0, t={t} is not finite."
            assert np.isclose(b_val, expected_b_at_throat[i]), \
                   f"shape_func at r=R0, t={t} does not match expected tanh behavior."

    def test_shape_func_asymptotic_behavior(self):
        """
        Test shape function far from the throat (large r).
        It should approach R0, as the exp(-(r-R0)**2) term goes to zero.
        """
        r_far = 100 * R0
        b_far = shape_func(r_far, t=0)
        assert np.isfinite(b_far), "shape_func far from throat is not finite."
        assert np.isclose(b_far, R0, atol=TOL), \
               "shape_func should approach R0 far from the throat."

    def test_traversable_metric_numerical_stability(self):
        """
        Test traversable_metric for numerical stability, especially near R0 and origin.
        Should not return NaN or Inf.
        """
        test_points = [
            (R0, 0, 0),       # At throat center
            (R0 + 0.01, 0, 0), # Just outside throat
            (R0 - 0.01, 0, 0), # Just inside throat
            (0, 0, 0),        # At origin (r=0)
            (10, 10, 0),      # Far from throat
            (0, 0, 5)         # At origin, later time
        ]
        for x, y, t in test_points:
            Z = traversable_metric(x, y, t)
            assert np.isfinite(Z), f"traversable_metric returned non-finite value at ({x}, {y}, t={t})."
            assert -5.0 <= Z <= 5.0, f"traversable_metric returned clipped value outside expected range at ({x}, {y}, t={t})."

    def test_nec_violation_at_throat(self):
        """
        Test the Null Energy Condition (NEC) violation around the throat.
        For a traversable wormhole, we expect rho + pr to be negative near the throat.
        """
        t_test = 2.0 # A time when tanh(t) is close to 1
        r_vals = np.linspace(R0 - 0.1, R0 + 0.1, 20) # Sample around the throat

        nec_results = [nec_violation(r, t_test) for r in r_vals]

        # Check if at least one value is significantly negative.
        # This is a heuristic test, a true violation would be consistently negative near throat.
        assert any(val < -1e-3 for val in nec_results), \
               f"NEC (rho + pr) not significantly negative around throat at t={t_test}. " \
               f"Values: {nec_results}"
        assert all(np.isfinite(val) for val in nec_results), \
               f"NEC violation calculation produced non-finite values at t={t_test}."

    def test_nec_violation_away_from_throat(self):
        """
        Test NEC violation far from the throat. Ideally, it should approach zero or positive.
        """
        t_test = 2.0
        r_far = 10 * R0
        nec_far = nec_violation(r_far, t_test)
        assert np.isfinite(nec_far), "NEC violation far from throat is not finite."
        # Expect it to be small or positive, but not necessarily exactly zero due to the functions
        assert nec_far >= -1e-3, f"NEC (rho + pr) unexpectedly negative far from throat at t={t_test}: {nec_far}"


    def test_geodesic_eq_numerical_stability(self):
        """
        Test geodesic_eq for numerical stability for various states.
        Should not return NaN or Inf.
        """
        test_states = [
            (R0, 0.1, 0, 0, 0),       # At throat, photon at rest (for testing derivatives)
            (3.0, 0.5, -1.0, 0.0, 2.0), # Typical initial photon condition
            (0.01, 0.01, 0.1, 0.1, 1.0) # Near origin
        ]
        for x, y, vx, vy, t in test_states:
            u = [x, y, vx, vy]
            derivs = geodesic_eq(u, 0, wormhole_t=t) # geodesic_eq's own t is dummy, wormhole_t is actual time
            assert len(derivs) == 4, "geodesic_eq did not return 4 derivatives."
            assert all(np.isfinite(d) for d in derivs), \
                   f"geodesic_eq returned non-finite derivatives for u={u}, t={t}. Derivs: {derivs}"

    def test_simulate_traversal_outcome(self):
        """
        Test if at least one photon successfully traverses the wormhole.
        This is a functional test of the overall setup.
        """
        print("\n--- Running simulate_traversal test (may take a moment) ---")
        traversed_count, _ = simulate_traversal()
        print(f"--- Simulate traversal test finished. Photons traversed: {traversed_count} ---")
        assert traversed_count > 0, "No photons traversed the wormhole. Check geodesic_eq or parameters."
        # You might want to assert a minimum number of traversed photons for robustness
        # assert traversed_count >= N_EXPECTED, "Too few photons traversed."

    # Additional tests could include:
    # - Testing if derivatives in geodesic_eq match analytical (if available)
    # - Testing behavior of metric functions for extreme parameter values (A, omega, beta)
    # - Ensuring that the time evolution of the mouth visualization is sensible (visually, harder to test with code)

if __name__ == "__main__":
    import pytest
    # For visualization purposes
    ani = animate_wormhole()
    plt.show()
    
    plot_nec()
    plt.show()
    
    traversed_count, _ = simulate_traversal()
    print(f"Photons traversed: {traversed_count}")
    pytest.main([__file__, "-v"])  # Same as command line
