import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, pandas as pd, numpy as np
d = pd.read_csv("pairs.csv")
x = d["dS_2day_cfs"].values; y = d["reg_minus_unreg_1hr"].values
fig, ax = plt.subplots(figsize=(7.2, 5.2), dpi=200)
ax.scatter(x, y, s=34, color="#2b6cb0", zorder=3)
for wy, xi, yi in zip(d["WY"], x, y):
    ax.annotate(str(int(wy)), (xi, yi), textcoords="offset points",
                xytext=(5, 5), fontsize=7.5, color="#333333")
xs = np.linspace(x.min(), x.max(), 50)
ax.plot(xs, -0.869 * xs - 3836, color="#c53030", lw=1.8, zorder=2,
        label="REG - UNREG = -0.869 x (2-day storage change) - 3,836")
ax.set_xlabel("Maximum 2-day Mossyrock storage change (mean cfs over window)", fontsize=10)
ax.set_ylabel("Regulated minus unregulated peak (cfs)", fontsize=10)
ax.set_title("R\u00b2 = 0.871    Standard error = 5,528 cfs    n = 17", fontsize=10)
ax.grid(alpha=0.25, lw=0.6)
ax.legend(fontsize=8.5, loc="lower left", frameon=True)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: format(int(v), ",")))
ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: format(int(v), ",")))
fig.tight_layout(); fig.savefig("reg2day.png", dpi=200, facecolor="white")
