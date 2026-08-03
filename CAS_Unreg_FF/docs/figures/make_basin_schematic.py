import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrow

fig, ax = plt.subplots(figsize=(9.0, 4.6), dpi=200)
ax.set_xlim(-2, 102); ax.set_ylim(0, 55); ax.axis("off")

blue = "#2b6cb0"; grey = "#4a5568"; red = "#c53030"

# main stem Cowlitz: headwaters (right) to Columbia (left)
main_x = [92, 84, 76, 68, 60, 52, 44, 36, 28, 20, 12, 6]
main_y = [40, 38, 36, 34, 33, 31, 29, 28, 26, 25, 24, 23]
ax.plot(main_x, main_y, color=blue, lw=3.0, solid_capstyle="round", zorder=2)

# Columbia River
ax.plot([6, 6], [46, 4], color=blue, lw=5.0, solid_capstyle="round", zorder=1)
ax.text(10.5, 8, "Columbia\nRiver", color=blue, fontsize=9, ha="center", style="italic")

# Toutle River tributary
ax.plot([22, 26, 32, 40], [8, 13, 18, 24], color=blue, lw=2.0, zorder=2)
ax.text(23.5, 6.0, "Toutle River", color=blue, fontsize=8.5, style="italic")

# Cispus tributary
ax.plot([80, 84, 90], [22, 28, 33], color=blue, lw=1.8, zorder=2)
ax.text(83.0, 19.5, "Cispus River", color=blue, fontsize=8.5, style="italic")

def dam(x, y, label, sub, color=grey):
    ax.plot([x, x], [y - 3.2, y + 3.2], color=color, lw=3.2, zorder=4)
    ax.text(x, y + 4.6, label, ha="center", fontsize=9.5, fontweight="bold")
    ax.text(x, y + 10.2, sub, ha="center", fontsize=8.0, color="#555555")
    ax.annotate("", xy=(x, y + 4.2), xytext=(x, y + 9.6),
                arrowprops=dict(arrowstyle="-", color="#bbbbbb", lw=0.8))

# reservoirs as pools
def pool(x0, x1, ymid, h, label):
    poly = Polygon([(x0, ymid - h), (x1, ymid - h * 0.55),
                    (x1, ymid + h * 0.55), (x0, ymid + h)],
                   closed=True, facecolor="#bee3f8", edgecolor=blue, lw=0.8, zorder=1)
    ax.add_patch(poly)
    ax.text((x0 + x1) / 2.0, ymid - h - 3.0, label, ha="center",
            fontsize=8.0, color=blue, style="italic")

pool(64, 78, 34.5, 3.6, "Riffe Lake")
pool(53.0, 59.0, 31.6, 2.2, "Mayfield Lake")

dam(64, 34.5, "Mossyrock Dam", "completed 1968; primary\nflood storage on the basin")
dam(52.5, 31.5, "Mayfield Dam", "completed 1963\nre-regulating, minor storage")
dam(84, 37.5, "Cowlitz Falls Dam", "1994; run-of-river,\nnegligible storage")

# gages
def gage(x, y, label, sub, dy=-4.0, dy2=-7.4):
    ax.plot([x], [y], marker="v", color=red, markersize=10, zorder=5)
    ax.text(x, y + dy, label, ha="center", fontsize=9.0, fontweight="bold", color=red)
    ax.text(x, y + dy2, sub, ha="center", fontsize=7.8, color="#555555")

gage(46.0, 29.6, "Mayfield gage", "USGS 14238000", dy=-8.6, dy2=-11.8)
gage(16.0, 24.6, "Castle Rock", "USGS 14243000")

ax.annotate("", xy=(4.5, 26), xytext=(11, 24),
            arrowprops=dict(arrowstyle="->", color=blue, lw=1.4))
ax.text(96, 43, "headwaters\n(Mount Rainier /\nGoat Rocks)", ha="center",
        fontsize=8.0, color="#555555")
ax.text(50, 51.0, "Cowlitz River Basin \u2014 schematic (not to scale); flow is right to left",
        ha="center", fontsize=10.5, fontweight="bold")

fig.tight_layout()
fig.savefig("basin_schematic.png", dpi=200, bbox_inches="tight",
            facecolor="white")
