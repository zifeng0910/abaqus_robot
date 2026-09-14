"""Create the three Python-only visual identity figures for the frozen CAD."""

from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
STEM = "Robot_L2300_D0815_WallWobble"
BLUE = "#2878B5"
GREEN = "#3A923A"
ORANGE = "#D9822B"
INK = "#252525"
LIGHT = "#D9D9D9"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "savefig.dpi": 400,
})


def load():
    profile = pd.read_csv(HERE/(STEM+"_profile.csv"))
    geometry = json.loads((HERE/(STEM+"_geometry.json")).read_text())
    props = json.loads((HERE/(STEM+"_mass_properties.json")).read_text())
    return profile, geometry, props


def save_bundle(fig, stem):
    fig.savefig(HERE/(stem+".png"), dpi=400, bbox_inches="tight")
    fig.savefig(HERE/(stem+".svg"), bbox_inches="tight")
    fig.savefig(HERE/(stem+".pdf"), bbox_inches="tight")
    fig.savefig(HERE/(stem+".tiff"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})


def radius_at(profile, x, length):
    head = profile[profile.region == "HEAD"]
    xp = head.x_mm.to_numpy()
    fp = head.radius_mm.to_numpy()
    if not np.all(np.diff(xp) >= 0):
        raise ValueError("HEAD axial profile must be monotone for interpolation")
    index = np.clip(np.searchsorted(xp, x, side="right"), 1, len(xp)-1)
    fraction = (x-xp[index-1])/np.maximum(xp[index]-xp[index-1], 1e-15)
    head_radius = fp[index-1] + fraction*(fp[index]-fp[index-1])
    return np.where(x <= head.x_mm.max(),
                    head_radius,
                    np.where(x <= length, 0.4075, np.nan))


def dimension(ax, start, end, y, text, color=INK):
    ax.annotate("", xy=(start, y), xytext=(end, y),
                arrowprops={"arrowstyle": "<->", "color": color, "lw": .9})
    ax.text((start+end)/2, y+0.035, text, ha="center", va="bottom", color=color)


def plot_profile(profile, geometry, props):
    head = profile[profile.region == "HEAD"]
    xh = geometry["head"]["axial_extent_mm"]
    length = geometry["L_total_mm"]
    radius = geometry["R_body_mm"]
    com = props["center_of_mass"]["value"][0]
    fig, ax = plt.subplots(figsize=(7.2, 3.1))
    ax.fill_between(head.x_mm, -head.radius_mm, head.radius_mm, color=BLUE, alpha=.25)
    ax.plot(head.x_mm, head.radius_mm, color=BLUE, lw=1.8)
    ax.plot(head.x_mm, -head.radius_mm, color=BLUE, lw=1.8)
    ax.fill_between([xh, length], [-radius, -radius], [radius, radius], color=GREEN, alpha=.20)
    ax.plot([xh, length], [radius, radius], color=GREEN, lw=1.8)
    ax.plot([xh, length], [-radius, -radius], color=GREEN, lw=1.8)
    ax.plot([length, length], [-radius, radius], color=ORANGE, lw=1.8)
    ax.scatter([com], [0], marker="x", s=35, linewidth=1.4, color=INK, zorder=5)
    ax.text(.085, .18, "HEAD\nBezier-5", color=BLUE, ha="center")
    ax.text((xh+length)/2, .08, "straight cylinder", color=GREEN, ha="center")
    ax.text(length-.02, -.12, "flat TAIL", color=ORANGE, ha="right")
    ax.text(com, -.10, "COM", color=INK, ha="center")
    dimension(ax, 0, length, -.61, "L = 2.300000 mm")
    ax.annotate("", xy=(0, .53), xytext=(xh, .53),
                arrowprops={"arrowstyle": "<->", "color": BLUE, "lw": .9})
    ax.text(.34, .56, "HEAD = 0.261495 mm", ha="left", va="bottom", color=BLUE)
    ax.annotate("D = 0.815000 mm", xy=(1.75, radius), xytext=(1.75, -.34),
                ha="center", va="center", rotation=90,
                rotation_mode="anchor",
                arrowprops={"arrowstyle": "<->", "color": GREEN, "lw": .9}, color=GREEN)
    ax.set(xlabel="Global X (mm)", ylabel="Radial coordinate (mm)",
           title="Authoritative L2300 D0815 FreeCAD profile")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-.08, 2.39)
    ax.set_ylim(-.69, .69)
    ax.grid(color=LIGHT, lw=.5, alpha=.65)
    fig.tight_layout()
    save_bundle(fig, STEM+"_profile")
    plt.close(fig)


def plot_preview(profile, geometry, props):
    length = geometry["L_total_mm"]
    x = np.linspace(0, length, 220)
    radius = radius_at(profile, x, length)
    phi = np.linspace(0, 2*np.pi, 72)
    X, P = np.meshgrid(x, phi)
    R = np.tile(radius, (len(phi), 1))
    Y, Z = R*np.cos(P), R*np.sin(P)
    fig = plt.figure(figsize=(7.2, 4.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z, color="#86A9C2", linewidth=0, antialiased=True, shade=True)
    com = props["center_of_mass"]["value"][0]
    ax.scatter([com], [0], [0], color=ORANGE, s=24, depthshade=False)
    fig.suptitle("FreeCAD BRep identity preview", y=.97)
    fig.text(.5, .895,
             "HEAD: Bezier-5, 0.261495 mm   |   BODY: straight, 2.038505 mm   |   TAIL: flat",
             ha="center", color=INK)
    fig.text(.5, .855, "L = 2.300000 mm   |   D = 0.815000 mm   |   COM X = 1.206019 mm",
             ha="center", color=ORANGE)
    ax.set(xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)")
    ax.set_box_aspect((2.3, .815, .815))
    ax.view_init(elev=22, azim=-58)
    ax.grid(False)
    fig.tight_layout(rect=(0, 0, 1, .80))
    save_bundle(fig, STEM+"_preview")
    plt.close(fig)


def plot_parameterization(profile, geometry):
    head = profile[profile.region == "HEAD"]
    xh = geometry["head"]["axial_extent_mm"]
    radius = geometry["R_body_mm"]
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    for length, color, label, yoff in (
        (1.8, BLUE, "L1.800", 0.0),
        (2.3, ORANGE, "L2.300", 0.0),
    ):
        ax.plot(head.x_mm, head.radius_mm+yoff, color=color, lw=1.6)
        ax.plot([xh, length], [radius+yoff, radius+yoff], color=color, lw=1.6, label=label)
        ax.plot([length, length], [-radius+yoff, radius+yoff], color=color, lw=1.2)
    ax.axvspan(1.8, 2.3, color=GREEN, alpha=.16)
    ax.annotate("+0.500 mm straight cylinder", xy=(2.05, .4075), xytext=(2.05, .62),
                ha="center", arrowprops={"arrowstyle": "->", "color": GREEN}, color=GREEN)
    ax.annotate("identical HEAD", xy=(.12, .28), xytext=(.45, .12),
                ha="center", color=INK,
                arrowprops={"arrowstyle": "->", "color": INK, "lw": .8})
    ax.set(xlabel="Global X (mm)", ylabel="Upper profile radius (mm)",
           title="Length parameterization changes only the straight body")
    ax.set_xlim(-.05, 2.38)
    ax.set_ylim(-.48, .72)
    ax.grid(color=LIGHT, lw=.5, alpha=.65)
    ax.legend(loc="lower left")
    fig.tight_layout()
    save_bundle(fig, "L1800_vs_L2300_parameterization")
    plt.close(fig)


def main():
    profile, geometry, props = load()
    plot_profile(profile, geometry, props)
    plot_preview(profile, geometry, props)
    plot_parameterization(profile, geometry)


if __name__ == "__main__":
    main()
