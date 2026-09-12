"""Publication-style one-variable comparison of zeta=0.20 and zeta=0.50."""
from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REF = '#5B6573'
CAND = '#D55E00'
ACCENT = '#0072B2'
GRID = '#D9DDE2'

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


def save(fig, stem):
    fig.savefig(HERE / f'{stem}.png', dpi=600, bbox_inches='tight')
    fig.savefig(HERE / f'{stem}.tiff', dpi=600, bbox_inches='tight')
    fig.savefig(HERE / f'{stem}.pdf', bbox_inches='tight')
    fig.savefig(HERE / f'{stem}.svg', bbox_inches='tight')
    plt.close(fig)


def impact_lines(ax, ref, cand, column, ylabel, zero=True):
    for frame, color, label in [(ref, REF, 'ζ=0.20'), (cand, CAND, 'ζ=0.50')]:
        mask = frame.time_s.between(.001000, .001010)
        ax.plot(frame.loc[mask, 'time_s'] * 1e3, frame.loc[mask, column],
                color=color, lw=1.25, label=label)
    if zero:
        ax.axhline(0, color=GRID, lw=.7, zorder=0)
    ax.set(xlabel='Time (ms)', ylabel=ylabel, xlim=(1.000, 1.010))
    ax.legend(loc='best')


def main():
    sr = json.loads((HERE / 'normaldamp020_summary.json').read_text())
    sc = json.loads((HERE / 'normaldamp050_summary.json').read_text())
    hr = pd.read_csv(HERE / 'normaldamp020_first_impact_history.csv')
    hc = pd.read_csv(HERE / 'normaldamp050_first_impact_history.csv')
    fr = pd.read_csv(HERE / 'normaldamp020_contact_force_history.csv')
    fc = pd.read_csv(HERE / 'normaldamp050_contact_force_history.csv')
    er = pd.read_csv(HERE / 'normaldamp020_energy_history.csv')
    ec = pd.read_csv(HERE / 'normaldamp050_energy_history.csv')

    fig, ax = plt.subplots(figsize=(3.5039, 2.45), constrained_layout=True)
    vals = [sr['e_n'], sc['e_n']]
    ax.axhspan(.30, .65, color='#009E73', alpha=.11, label='Provisional window')
    ax.plot([.20, .50], vals, color=ACCENT, lw=1.2, marker='o', ms=4)
    for x, y in zip([.20, .50], vals):
        ax.text(x, y + .025, f'{y:.3f}', ha='center', va='bottom')
    ax.set(xlabel='Normal critical damping fraction, ζ', ylabel='Normal restitution, e_n',
           xlim=(.14, .56), ylim=(.25, .90), xticks=[.20, .50])
    ax.legend(loc='lower left')
    save(fig, 'normaldamp020_vs_050_restitution')

    fig, ax = plt.subplots(figsize=(3.5039, 2.30), constrained_layout=True)
    impact_lines(ax, hr, hc, 'gap_um', 'Minimum node-to-wall gap (µm)')
    ax.fill_between([1.000, 1.010], -1, 0, color=GRID, alpha=.25, zorder=0)
    save(fig, 'normaldamp020_vs_050_gap')

    fig, ax = plt.subplots(figsize=(3.5039, 2.30), constrained_layout=True)
    impact_lines(ax, fr, fc, 'Fwall_norm_N', 'Wall-force magnitude (N)', zero=False)
    ax.set_ylim(bottom=0)
    save(fig, 'normaldamp020_vs_050_wall_force')

    fig, ax = plt.subplots(figsize=(3.5039, 2.45), constrained_layout=True)
    for frame, summary, color, label in [(er, sr, REF, 'ζ=0.20'), (ec, sc, CAND, 'ζ=0.50')]:
        mask = frame.time_s.between(.001000, .001010)
        anchor = int(np.argmin(np.abs(frame.time_s.to_numpy() - summary['first_contact_s'])))
        delta = (frame.Ktotal_J - frame.Ktotal_J.iloc[anchor]) * 1e6
        ax.plot(frame.loc[mask, 'time_s'] * 1e3, delta[mask], color=color, lw=1.25, label=label)
    ax.axhline(0, color=GRID, lw=.7)
    ax.set(xlabel='Time (ms)', ylabel='Change in total kinetic energy (µJ)', xlim=(1.000, 1.010))
    ax.legend(loc='best')
    save(fig, 'normaldamp020_vs_050_energy')

    fig, ax = plt.subplots(figsize=(3.5039, 2.30), constrained_layout=True)
    recoil = [sr['deltaV_COM_norm_mm_s'], sc['deltaV_COM_norm_mm_s']]
    bars = ax.bar(['ζ=0.20', 'ζ=0.50'], recoil, color=[REF, CAND], width=.56)
    ax.bar_label(bars, labels=[f'{value:.1f}' for value in recoil], padding=3, fontsize=7)
    ax.set(ylabel='First-impact COM recoil (mm s^-1)', ylim=(0, max(recoil) * 1.18))
    ax.text(.98, .97,
            f'Reduction: {sc["COM_recoil_reduction_percent"]:.1f}% (secondary diagnostic)',
            transform=ax.transAxes, ha='right', va='top')
    save(fig, 'normaldamp020_vs_050_COM_recoil')


if __name__ == '__main__':
    main()
