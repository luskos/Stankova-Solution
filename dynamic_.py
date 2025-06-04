import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# =============================================
# 1. Traversable Wormhole Metric
# =============================================
R0 = 1.0       # Throat radius
A = 0.5        # Flare control
omega = 0.3    # Oscillation frequency
beta = 0.1     # Asymmetry parameter

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
# 2. Dynamic Visualization (Mouth Opening)
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
    plt.show()
    return ani

# =============================================
# 3. Geodesics for Traversing Photons
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
    plt.show()


def simulate_traversal():
    plt.figure(figsize=(8, 8))
    t_span = np.linspace(0, 50, 3000)
    
    traversed_count = 0

    # Initial conditions (photon aimed at throat)
    for impact_param in np.linspace(-0.5, 0.5, 15):  # more photons for better sampling
        sol = odeint(geodesic_eq, [3.0, impact_param, -1.0, 0.0], t_span, args=(2.0,))
        plt.plot(sol[:, 0], sol[:, 1], label=f"b={impact_param:.2f}")
        
        # Check if photon crossed x=0
        crossing_indices = np.where(np.diff(np.sign(sol[:, 0])))[0]
        if len(crossing_indices) > 0:
            print(f"✅ Photon with b = {impact_param:.2f} traversed.")
            traversed_count += 1

    plt.gca().add_patch(plt.Circle((0, 0), R0, color='red', alpha=0.2, label='Throat'))
    plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    plt.title("Photon Traversal Through Wormhole")
    plt.xlabel("X"), plt.ylabel("Y")
    plt.axis('equal'), plt.grid()
    plt.legend()
    plt.show()

    print(f"\n🔎 Total traversing photons: {traversed_count}")


# =============================================
# Run Tests
# =============================================
animate_wormhole()  # Shows mouth opening/closing
simulate_traversal()  # Photons passing through throat
plot_nec(t=2.0)  # Или друго време, което ти е интересно