import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

from matplotlib.colors import Normalize
from matplotlib.path import Path
from scipy.interpolate import griddata


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
#
# x, y, u, v, is_boundary
#
# is_boundary = 1 -> boundary edge
# is_boundary = 0 -> internal edge
# ============================================================

edge = np.genfromtxt(
    "edge_data.csv",
    delimiter=",",
    names=True
)

xe = edge["x"]
ye = edge["y"]

u = edge["u"]
v = edge["v"]

is_boundary = edge["is_boundary"].astype(bool)


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
# Find boundary edges from the TRIANGULAR MESH
#
# We reconstruct every triangle edge.
#
# Each edge is represented using the coordinates of its
# two endpoints.
# ============================================================

all_edges = {}

for i in range(ncell):

    vertices = [
        (x1[i], y1[i]),
        (x2[i], y2[i]),
        (x3[i], y3[i])
    ]

    triangle_edges = [
        (vertices[0], vertices[1]),
        (vertices[1], vertices[2]),
        (vertices[2], vertices[0])
    ]

    for p1, p2 in triangle_edges:

        # Sort endpoints so that the same edge has the
        # same key regardless of orientation.
        key = tuple(sorted([
            tuple(p1),
            tuple(p2)
        ]))

        if key not in all_edges:
            all_edges[key] = 0

        all_edges[key] += 1


# ============================================================
# A mesh boundary edge occurs only once.
#
# Interior mesh edges occur twice.
# ============================================================

boundary_edges = []

for edge_key, count in all_edges.items():

    if count == 1:

        p1 = np.array(edge_key[0])
        p2 = np.array(edge_key[1])

        boundary_edges.append(
            (p1, p2)
        )


print("\nMesh boundary information")
print("-------------------------")

print(
    f"Number of boundary edges found: "
    f"{len(boundary_edges)}"
)


# ============================================================
# Build connected boundary loops
#
# This is important because there can be:
#
#   1. Outer CFD boundary
#   2. Airfoil/object boundary
#
# We find all closed boundary loops automatically.
# ============================================================

def build_boundary_loops(boundary_edges):

    # Convert coordinates to tuples for dictionary lookup

    adjacency = {}

    for p1, p2 in boundary_edges:

        p1 = tuple(p1)
        p2 = tuple(p2)

        if p1 not in adjacency:
            adjacency[p1] = []

        if p2 not in adjacency:
            adjacency[p2] = []

        adjacency[p1].append(p2)
        adjacency[p2].append(p1)

    unused = set()

    for p1, p2 in boundary_edges:

        unused.add(
            tuple(sorted([
                tuple(p1),
                tuple(p2)
            ]))
        )

    loops = []

    while unused:

        edge0 = next(iter(unused))

        start = edge0[0]
        current = start

        loop = [start]

        previous = None

        while True:

            candidates = adjacency[current]

            next_point = None

            for candidate in candidates:

                edge_key = tuple(sorted([
                    current,
                    candidate
                ]))

                if edge_key in unused:

                    if candidate != previous or len(candidates) == 1:
                        next_point = candidate
                        break

            if next_point is None:
                break

            edge_key = tuple(sorted([
                current,
                next_point
            ]))

            unused.remove(edge_key)

            previous = current
            current = next_point

            if current == start:
                break

            loop.append(current)

        if len(loop) >= 3:

            loops.append(
                np.array(loop)
            )

    return loops


boundary_loops = build_boundary_loops(
    boundary_edges
)


print(
    f"Number of closed boundary loops: "
    f"{len(boundary_loops)}"
)


# ============================================================
# Identify the OBJECT / AIRFOIL LOOP
#
# Usually:
#
#   largest loop  -> outer domain
#   smaller loop  -> airfoil
#
# We therefore calculate the area of every loop.
# ============================================================

def polygon_area(points):

    x = points[:, 0]
    y = points[:, 1]

    return 0.5 * abs(
        np.sum(
            x * np.roll(y, -1)
            -
            y * np.roll(x, -1)
        )
    )


loop_areas = []

for loop in boundary_loops:

    loop_areas.append(
        polygon_area(loop)
    )


print("\nBoundary loop areas")
print("-------------------")

for i, area in enumerate(loop_areas):

    print(
        f"Loop {i}: area = {area:.6g}"
    )


# ============================================================
# Assume the largest loop is the OUTER CFD boundary.
#
# The smaller closed loop is the airfoil/object.
# ============================================================

if len(boundary_loops) < 2:

    raise RuntimeError(
        "\nCould not find an internal object boundary.\n"
        "Make sure the airfoil is represented by a closed "
        "boundary in the CFD mesh."
    )


outer_index = np.argmax(loop_areas)

object_indices = [
    i for i in range(len(boundary_loops))
    if i != outer_index
]


# ============================================================
# If there is only one object, use it.
#
# If there are multiple internal objects, all of them
# will be masked.
# ============================================================

object_loops = [
    boundary_loops[i]
    for i in object_indices
]


print(
    f"\nObject boundaries detected: "
    f"{len(object_loops)}"
)


# ============================================================
# Create the plotting grid
# ============================================================

xmin = points_x.min()
xmax = points_x.max()

ymin = points_y.min()
ymax = points_y.max()

Nx = 1000
Ny = 1000

gx = np.linspace(
    xmin,
    xmax,
    Nx
)

gy = np.linspace(
    ymin,
    ymax,
    Ny
)

GX, GY = np.meshgrid(
    gx,
    gy
)


# ============================================================
# Interpolate velocity
# ============================================================

GU = griddata(
    (xe, ye),
    u,
    (GX, GY),
    method="linear"
)

GV = griddata(
    (xe, ye),
    v,
    (GX, GY),
    method="linear"
)


# ============================================================
# Velocity magnitude
# ============================================================

GSPEED = np.sqrt(
    GU**2 +
    GV**2
)


# ============================================================
# Interpolate pressure from cell centers
# ============================================================

xc = (
    x1 +
    x2 +
    x3
) / 3.0

yc = (
    y1 +
    y2 +
    y3
) / 3.0


GP = griddata(
    (xc, yc),
    P,
    (GX, GY),
    method="linear"
)


# ============================================================
# Build OBJECT MASK
#
# This is the important part.
#
# We are NOT drawing an airfoil on top.
#
# We identify the actual closed boundary from the CFD mesh
# and determine which plotting-grid points lie inside it.
# ============================================================

inside_object = np.zeros(
    GX.shape,
    dtype=bool
)


grid_points = np.column_stack([
    GX.ravel(),
    GY.ravel()
])


for loop in object_loops:

    object_path = Path(loop)

    inside = object_path.contains_points(
        grid_points
    )

    inside = inside.reshape(
        GX.shape
    )

    inside_object |= inside


# ============================================================
# Mask the OBJECT from pressure
# ============================================================

GP_plot = np.where(
    inside_object,
    np.nan,
    GP
)


# ============================================================
# Mask invalid velocity data
# ============================================================

valid_velocity = (
    np.isfinite(GU) &
    np.isfinite(GV)
)


GU_plot = np.where(
    valid_velocity & (~inside_object),
    GU,
    np.nan
)

GV_plot = np.where(
    valid_velocity & (~inside_object),
    GV,
    np.nan
)

GSPEED_plot = np.where(
    valid_velocity & (~inside_object),
    GSPEED,
    np.nan
)


# ============================================================
# Figure
# ============================================================

fig, ax = plt.subplots(
    figsize=(12, 8)
)


# ============================================================
# Smooth pressure contours
# ============================================================

pressure_plot = ax.contourf(
    GX,
    GY,
    GP_plot,
    levels=100,
    cmap="viridis"
)


# ============================================================
# Pressure colorbar
# ============================================================

cbar1 = fig.colorbar(
    pressure_plot,
    ax=ax,
    pad=0.02
)

cbar1.set_label(
    "Pressure"
)


# ============================================================
# Streamlines
#
# The object region contains NaN.
#
# Therefore streamplot cannot draw streamlines through
# the airfoil.
# ============================================================

stream = ax.streamplot(
    gx,
    gy,
    GU_plot,
    GV_plot,

    color=GSPEED_plot,

    cmap="plasma",

    density=3.0,

    linewidth=1.2,

    arrowsize=1.0
)


# ============================================================
# Velocity colorbar
# ============================================================

cbar2 = fig.colorbar(
    stream.lines,
    ax=ax,
    pad=0.08
)

cbar2.set_label(
    "Velocity Magnitude (m/s)"
)


# ============================================================
# Plot the ACTUAL CFD boundary
#
# This is only the boundary extracted from the mesh.
# We are NOT generating an artificial NACA airfoil.
# ============================================================

for loop in object_loops:

    loop_closed = np.vstack([
        loop,
        loop[0]
    ])

    ax.plot(
        loop_closed[:, 0],
        loop_closed[:, 1],

        linewidth=2.0,

        color="black",

        zorder=20
    )


# ============================================================
# Formatting
# ============================================================

ax.set_xlabel("x")

ax.set_ylabel("y")

ax.set_title(
    "CFD Flow Field: Pressure Contours + Streamlines"
)

ax.set_aspect(
    "equal"
)

plt.tight_layout()

plt.show()