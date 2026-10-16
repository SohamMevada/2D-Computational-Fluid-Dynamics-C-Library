import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

filename = "edges_after_cleanup.txt"

fig, ax = plt.subplots(figsize=(10, 8))

with open(filename, "r") as f:

    for line in f:

        data = line.split()

        if len(data) != 8:
            continue

        edge_id = int(data[0])

        x1 = float(data[1])
        y1 = float(data[2])
        x2 = float(data[3])
        y2 = float(data[4])

        isBounEdge = int(data[5])
        isPdef = int(data[6])
        isVdef = int(data[7])

        # --------------------------------------------------
        # Determine boundary type
        # --------------------------------------------------

        if not isBounEdge:
            boundary_type = "internal"
            linestyle = "-"
        elif isPdef and isVdef:
            boundary_type = "both"
            linestyle = "-"
        elif isPdef:
            boundary_type = "pDef"
            linestyle = "-"
        elif isVdef:
            boundary_type = "vDef"
            linestyle = "-"
        else:
            boundary_type = "boundary_none"
            linestyle = "--"

        # --------------------------------------------------
        # Choose color
        # --------------------------------------------------

        if boundary_type == "internal":
            color = "black"

        elif boundary_type == "pDef":
            color = "red"

        elif boundary_type == "vDef":
            color = "blue"

        elif boundary_type == "both":
            color = "green"

        else:
            color = "orange"

        # --------------------------------------------------
        # Plot edge
        # --------------------------------------------------

        ax.plot(
            [x1, x2],
            [y1, y2],
            color=color,
            linestyle=linestyle,
            linewidth=2
        )

        # --------------------------------------------------
        # Edge ID at midpoint
        # --------------------------------------------------

        xm = (x1 + x2) / 2
        ym = (y1 + y2) / 2

        # ax.text(
        #     xm,
        #     ym,
        #     str(edge_id),
        #     fontsize=7,
        #     ha="center",
        #     va="center"
        # )


# ----------------------------------------------------------
# Legend
# ----------------------------------------------------------

legend_elements = [
    Line2D([0], [0],
           color="black",
           lw=2,
           label="Internal"),

    Line2D([0], [0],
           color="red",
           lw=2,
           label="Pressure defined (pDef)"),

    Line2D([0], [0],
           color="blue",
           lw=2,
           label="Velocity defined (vDef)"),

    Line2D([0], [0],
           color="green",
           lw=2,
           label="Both pDef + vDef"),

    Line2D([0], [0],
           color="orange",
           lw=2,
           linestyle="--",
           label="Boundary, neither defined")
]

ax.legend(handles=legend_elements)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Mesh Boundary Conditions")

ax.set_aspect("equal")
ax.grid(True)

plt.show()