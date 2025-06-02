import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

# =============================================
# 1. Дефиниране на параметри и метрика
# =============================================
R0 = 1.0  # Радиус на "хоризонта"
A = 0.5   # Сила на Инь-Ян асиметрията

def yin_yang_metric(x, y):
    r = np.sqrt(x**2 + y**2)
    θ = np.arctan2(y, x)
    Z = (1 - R0/r) * (A * np.sin(2*θ) * np.exp(-(r-R0)**2))
    return Z

# =============================================
# 2. Визуализация на топологията (3D и 2D)
# =============================================
def plot_topology():
    x = np.linspace(-3, 3, 100)
    y = np.linspace(-3, 3, 100)
    X, Y = np.meshgrid(x, y)
    Z = yin_yang_metric(X, Y)

    # 3D графика
    fig = plt.figure(figsize=(14, 6))
    ax1 = fig.add_subplot(121, projection='3d')
    surf = ax1.plot_surface(X, Y, Z, cmap='coolwarm', edgecolor='none')
    ax1.set_title("3D изглед на Инь-Ян топологията")
    fig.colorbar(surf, ax=ax1, shrink=0.5, label="Потенциал")

    # 2D контурна графика
    ax2 = fig.add_subplot(122)
    cntr = ax2.contourf(X, Y, Z, levels=20, cmap='coolwarm')
    ax2.add_patch(plt.Circle((0, 0), R0, color='black', alpha=0.3, label='Хоризонт'))
    ax2.set_title("2D контурна карта")
    ax2.set_xlabel("X"), ax2.set_ylabel("Y")
    ax2.legend(), plt.colorbar(cntr, ax=ax2)
    plt.tight_layout()
    plt.show()

# =============================================
# 3. Симулация на фотони (геодезични)
# =============================================
def geodesic_eq(u, t):
    x, y, vx, vy = u
    r = np.sqrt(x**2 + y**2)
    F = (A * (R0/r)**3) * np.sin(2*np.arctan2(y, x))
    ax = F * y/r 
    ay = -F * x/r
    return [vx, vy, ax, ay]

def simulate_photons():
    plt.figure(figsize=(10, 10))
    colors = plt.cm.rainbow(np.linspace(0, 1, 12))
    
    for i, θ in enumerate(np.linspace(0, 2*np.pi, 12, endpoint=False)):
        x0, y0 = 2.5 * np.cos(θ), 2.5 * np.sin(θ)
        vx0, vy0 = -0.2 * np.sin(θ), 0.2 * np.cos(θ)
        t = np.linspace(0, 50, 1000)
        sol = odeint(geodesic_eq, [x0, y0, vx0, vy0], t)
        plt.plot(sol[:, 0], sol[:, 1], color=colors[i], alpha=0.7, label=f"{np.degrees(θ):.0f}°")
    
    plt.gca().add_patch(plt.Circle((0, 0), R0, color='black', alpha=0.3))
    plt.title("Траектории на фотони\n(Отблъскване от бяла дупка + вихри)")
    plt.xlabel("X"), plt.ylabel("Y"), plt.legend()
    plt.axis('equal'), plt.grid()
    plt.show()

# =============================================
# 4. Анализ на енергийното разпределение
# =============================================
def analyze_energy():
    x = np.linspace(-3, 3, 200)
    y = np.linspace(-3, 3, 200)
    X, Y = np.meshgrid(x, y)
    Z = yin_yang_metric(X, Y)
    energy = np.abs(Z.flatten())
    radial_dist = np.sqrt(X.flatten()**2 + Y.flatten()**2)

    plt.figure(figsize=(14, 5))
    
    # Хистограма
    plt.subplot(121)
    plt.hist(energy, bins=50, color='purple', alpha=0.7, edgecolor='white')
    plt.xlabel("Енергийна плътност"), plt.ylabel("Честота")
    plt.title("Разпределение на енергията")
    
    # Енергия vs радиус
    plt.subplot(122)
    plt.scatter(radial_dist, energy, c=radial_dist, cmap="viridis", alpha=0.5)
    plt.xlabel("Разстояние от центъра"), plt.ylabel("Енергия")
    plt.colorbar(label="Радиус")
    plt.title("Зависимост енергия-разстояние")
    plt.tight_layout()
    plt.show()

# =============================================
# Главна функция
# =============================================
if __name__ == "__main__":
    print("=== Симулация на Инь-Ян топология като бяла дупка ===")
    plot_topology()
    simulate_photons()
    analyze_energy()