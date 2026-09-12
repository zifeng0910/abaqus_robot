"""Publication-style baseline-versus-probe figures from public CSV sources."""
from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
BASE = '#5B6573'
CAND = '#D55E00'
CONTACT = '#0072B2'

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.titlesize': 8,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'axes.linewidth': .7,
    'legend.frameon': False,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
})


def save(fig, name):
    fig.savefig(HERE / (name + '.png'), dpi=600, bbox_inches='tight')
    fig.savefig(HERE / (name + '.tiff'), dpi=600, bbox_inches='tight')
    fig.savefig(HERE / (name + '.pdf'), bbox_inches='tight')
    fig.savefig(HERE / (name + '.svg'), bbox_inches='tight')
    plt.close(fig)


def lines(ax, xb, yb, xc, yc, ylabel):
    mb = (xb >= .001) & (xb <= .00101)
    mc = (xc >= .001) & (xc <= .00101)
    ax.plot(xb[mb] * 1e3, np.asarray(yb)[mb], color=BASE, lw=1.2, label='Baseline: 0.055, tangent 1')
    ax.plot(xc[mc] * 1e3, np.asarray(yc)[mc], color=CAND, lw=1.2, label='Probe: 0.20, tangent 0')
    ax.axhline(0, color='#A7ADB5', lw=.6, zorder=0)
    ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(1.000, 1.010))
    ax.legend(loc='best')


def main():
    base = pd.read_csv(HERE / 'current_contact_normal_tangential_work.csv')
    cand = pd.read_csv(HERE / 'normaldamp020_first_impact_history.csv')
    sb = json.loads((HERE / 'current_contact_restitution_summary.json').read_text())
    sc = json.loads((HERE / 'normaldamp020_summary.json').read_text())

    fig, ax = plt.subplots(figsize=(3.5039, 2.2441), constrained_layout=True)
    lines(ax, base.time_s, base.gap_um, cand.time_s, cand.gap_um, 'Minimum node-to-wall gap (µm)')
    ax.fill_between([1.000, 1.010], -1, 0, color='#D9D9D9', alpha=.25, zorder=0)
    save(fig, 'baseline_vs_normaldamp020_gap')

    fig, ax = plt.subplots(figsize=(3.5039, 2.2441), constrained_layout=True)
    lines(ax, base.time_s, base.vn_mm_s, cand.time_s, cand.vn_mm_s,
          'Contact-point normal velocity (mm s-1)')
    ax.text(.97, .12, f'e_n: {sb["e_n"]:.3f} to {sc["e_n"]:.3f}',
            transform=ax.transAxes, ha='right', va='bottom')
    save(fig, 'baseline_vs_normaldamp020_normal_velocity')

    bf = pd.read_csv(AUDIT / 'first_impact_contact_force_history.csv')
    cf = pd.read_csv(HERE / 'normaldamp020_contact_force_history.csv')
    fig, ax = plt.subplots(figsize=(3.5039, 2.2441), constrained_layout=True)
    lines(ax, bf.time_s, np.linalg.norm(bf[['Ftotal_robot_x_N', 'Ftotal_robot_y_N', 'Ftotal_robot_z_N']], axis=1),
          cf.time_s, cf.Fwall_norm_N, 'Wall-force magnitude (N)')
    save(fig, 'baseline_vs_normaldamp020_wall_force')

    be = pd.read_csv(AUDIT / 'first_impact_energy_timeline.csv')
    ce = pd.read_csv(HERE / 'normaldamp020_energy_history.csv')
    fig, ax = plt.subplots(figsize=(3.5039, 2.4409), constrained_layout=True)
    for frame, color, label in [(be, BASE, 'Baseline'), (ce, CAND, 'Probe')]:
        mask = (frame.time_s >= .001) & (frame.time_s <= .00101)
        anchor = int(np.argmin(np.abs(frame.time_s.to_numpy() - .0010035)))
        y = (frame.Ktotal_J - frame.Ktotal_J.iloc[anchor]) * 1e6
        ax.plot(frame.time_s[mask] * 1e3, y[mask], color=color, lw=1.2, label=label + ' ΔK total')
    ax.axhline(0, color='#A7ADB5', lw=.6)
    ax.set(xlabel='Time (ms)', ylabel='Change in total kinetic energy (µJ)', xlim=(1.000, 1.010))
    ax.legend(loc='best')
    save(fig, 'baseline_vs_normaldamp020_energy')

    fig, ax = plt.subplots(figsize=(3.5039, 2.2441), constrained_layout=True)
    vals = [sb['deltaV_COM_norm_mm_s'], sc['deltaV_COM_norm_mm_s']]
    bars = ax.bar(['Baseline', 'Probe'], vals, color=[BASE, CAND], width=.58)
    ax.bar_label(bars, labels=[f'{v:.1f}' for v in vals], padding=3, fontsize=7)
    ax.set(ylabel='First-impact COM recoil (mm s-1)', ylim=(0, max(vals) * 1.18))
    ax.text(.5, max(vals) * 1.08, f'Reduction: {sc["COM_recoil_reduction_percent"]:.1f}%', ha='center')
    save(fig, 'baseline_vs_normaldamp020_COM_recoil')


if __name__ == '__main__':
    main()
