import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import Normalize


# ============================================================
# Read cell data
# ============================================================

cell = np.genfromtxt(
    "cell_data.csv",
    delimiter=",",
    names=True
)

x1 = cell["x1"]
y1 = cell["y1"]

x2 = cell["x2"]
y2 = cell["y2"]

x3 = cell["x3"]
y3 = cell["y3"]

P = cell["pressure"]


# ============================================================
# Read edge data
#
# Expected columns:
# x1, y1, x2, y2, u, v
# ============================================================

edge = np.genfromtxt(
    "edge_data.csv",
    delimiter=",",
    names=True
)

ex1 = edge["x"]
ey1 = edge["y"]

ex2 = edge["x"]
ey2 = edge["y"]

u = edge["u"]
v = edge["v"]


# ============================================================
# Edge midpoints
# ============================================================

xe = 0.5 * (ex1 + ex2)
ye = 0.5 * (ey1 + ey2)


# ============================================================
# Velocity magnitude
# ============================================================

speed = np.sqrt(u**2 + v**2)


# ============================================================
# Build triangular mesh
# ============================================================

points_x = np.concatenate([x1, x2, x3])
points_y = np.concatenate([y1, y2, y3])

ncell = len(cell)

triangles = np.column_stack([
    np.arange(ncell),
    np.arange(ncell) + ncell,
    np.arange(ncell) + 2 * ncell
])

triang = tri.Triangulation(
    points_x,
    points_y,
    triangles
)


# ============================================================
# Velocity diagnostics
# ============================================================

print("\nEdge velocity statistics")
print("------------------------")
print(f"Min velocity  : {np.min(speed):.6g} m/s")
print(f"Max velocity  : {np.max(speed):.6g} m/s")
print(f"Mean velocity : {np.mean(speed):.6g} m/s")

imax = np.argmax(speed)

print("\nMaximum edge velocity location")
print("------------------------------")
print(f"x = {xe[imax]:.6f}")
print(f"y = {ye[imax]:.6f}")
print(f"u = {u[imax]:.6f}")
print(f"v = {v[imax]:.6f}")
print(f"|V| = {speed[imax]:.6f} m/s")


# ============================================================
# Figure
# ============================================================

fig, ax = plt.subplots(figsize=(10, 7))


# ============================================================
# Pressure field
# ============================================================

pressure_plot = ax.tripcolor(
    triang,
    np.repeat(P, 3),
    shading="flat",
    cmap="viridis"
)


# ============================================================
# Mesh overlay
# ============================================================

ax.triplot(
    triang,
    linewidth=0.35,
    alpha=0.5
)


# ============================================================
# Edge-centered velocity field
# ============================================================

u_dir = np.zeros_like(u)
v_dir = np.zeros_like(v)

nonzero = speed > 1e-12

u_dir[nonzero] = u[nonzero] / speed[nonzero]
v_dir[nonzero] = v[nonzero] / speed[nonzero]


q = ax.quiver(
    xe,
    ye,
    u_dir,
    v_dir,
    speed,
    cmap="plasma",
    norm=Normalize(
        vmin=0,
        vmax=np.max(speed)
    ),
    angles="xy",
    scale_units="xy",
    scale=7,
    width=0.0025,
    pivot="mid"
)


# ============================================================
# Colorbars
# ============================================================

cbar1 = fig.colorbar(
    pressure_plot,
    ax=ax,
    pad=0.02
)

cbar1.set_label("Pressure")


cbar2 = fig.colorbar(
    q,
    ax=ax,
    pad=0.08
)

cbar2.set_label("Velocity Magnitude (m/s)")


# ============================================================
# Formatting
# ============================================================

ax.set_xlabel("x")
ax.set_ylabel("y")

ax.set_title(
    "CFD Solution: Pressure Field and Edge-Centered Velocity Field"
)

ax.set_aspect("equal")

plt.tight_layout()
plt.show()
