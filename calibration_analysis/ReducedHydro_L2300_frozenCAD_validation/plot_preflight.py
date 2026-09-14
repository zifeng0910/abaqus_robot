"""Plot the frozen-CAD static window and solver-surface preflight evidence."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7, "axes.linewidth": 0.8, "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})
BLUE, RED, GOLD, GREY = "#3973A5", "#B7463A", "#D99B2B", "#707070"


def save(fig, stem):
    fig.savefig(HERE / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(HERE / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(HERE / f"{stem}.tiff", dpi=600, bbox_inches="tight")
    fig.savefig(HERE / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main():
    static = pd.read_csv(HERE / "L2300_static_clearance_map.csv")
    fractions = pd.read_csv(HERE / "L2300_static_wall_support_fraction_vs_tilt.csv")
    regression = pd.read_csv(HERE / "L2300_CAD_vs_solver_surface_regression.csv")

    pivot = static.pivot(index="tilt_deg", columns="azimuth_deg", values="min_gap_um")
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    image = ax.imshow(pivot.values, aspect="auto", origin="lower", cmap="RdBu_r",
                      vmin=-200, vmax=200, extent=(-5, 355, 13.5, 39.5))
    ax.contour(pivot.columns, pivot.index, pivot.values, levels=[0, 20],
               colors=["black", GOLD], linewidths=[0.8, 1.0])
    ax.set(xlabel="Azimuth (deg)", ylabel="Directed static tilt (deg)")
    bar = fig.colorbar(image, ax=ax, pad=0.02)
    bar.set_label("Exact minimum gap (um)")
    ax.text(4, 38.2, "Black: 0 um   Gold: 20 um", fontsize=6,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85})
    save(fig, "L2300_static_clearance_map")

    fig, ax = plt.subplots(figsize=(3.7, 2.8))
    ax.plot(fractions.tilt_deg, fractions.wall_support_fraction, "o-", color=BLUE,
            label="Wall support")
    ax.plot(fractions.tilt_deg, fractions.opposing_20um_fraction, "s-", color=RED,
            label="Opposing at 20 um")
    ax.set(xlabel="Directed static tilt (deg)", ylabel="Azimuth fraction", ylim=(-0.03, 1.05))
    ax.axvspan(28, 33, color=GOLD, alpha=0.12, lw=0)
    ax.legend(loc="upper left")
    save(fig, "L2300_static_wall_support_fraction_vs_tilt")

    fig, ax = plt.subplots(figsize=(3.7, 3.2))
    good = regression.contact_classification_match & regression.wall_support_classification_match
    ax.scatter(regression.loc[good, "min_gap_um"], regression.loc[good, "mesh_min_gap_um"],
               s=11, color=GREY, alpha=0.45, edgecolors="none", label="All classifications match")
    ax.scatter(regression.loc[~good, "min_gap_um"], regression.loc[~good, "mesh_min_gap_um"],
               s=28, color=RED, edgecolors="white", linewidths=0.4, label="Hard-gate mismatch", zorder=3)
    lo = min(regression.min_gap_um.min(), regression.mesh_min_gap_um.min())
    hi = max(regression.min_gap_um.max(), regression.mesh_min_gap_um.max())
    ax.plot([lo, hi], [lo, hi], color="black", lw=0.8)
    ax.axhline(0, color=RED, lw=0.6, ls="--"); ax.axvline(0, color=RED, lw=0.6, ls="--")
    ax.axhline(20, color=GOLD, lw=0.6, ls=":"); ax.axvline(20, color=GOLD, lw=0.6, ls=":")
    ax.set(xlabel="Exact CAD gap (um)", ylabel="Solver-surface gap (um)", xlim=(lo-5, hi+5), ylim=(lo-5, hi+5))
    ax.legend(loc="upper left", fontsize=6)
    save(fig, "L2300_CAD_vs_solver_surface")


if __name__ == "__main__":
    main()
