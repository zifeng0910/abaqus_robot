"""First-impact energy transfer figure from complete increment histories."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
    'svg.fonttype': 'none', 'pdf.fonttype': 42,
    'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'axes.spines.right': False, 'axes.spines.top': False,
    'axes.linewidth': 0.8, 'legend.frameon': False,
})

d = pd.read_csv(HERE / 'first_impact_energy_timeline.csv')
w = pd.read_csv(HERE / 'first_impact_energy_windows.csv')
event = w.loc[w.window == 'first_impact_shoulders'].iloc[0]
t = d.time_s.to_numpy() * 1e3
assert np.all(np.diff(t) > 0)
active = d.Fwall_norm_N.to_numpy() > 1e-8
ta, tb = t[active][0], t[active][-1]
baseline = d.iloc[0]

fig, axes = plt.subplots(1, 2, figsize=(7.2047, 2.9921), gridspec_kw={'width_ratios': [1.55, 1.]})
ax = axes[0]
ax.axvspan(ta, tb, color='#D55E00', alpha=.12, lw=0)
ax.plot(t, d.Ktrans_J * 1e6, color='#0072B2', lw=1.35, label='Translational KE')
ax.plot(t, d.Krot_J * 1e6, color='#6B7280', lw=1.25, label='Rotational KE')
ax.plot(t, d.Ktotal_J * 1e6, color='#1F2937', lw=1.0, ls='--', label='Total KE')
ax.set_xlim(.998, 1.012)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Kinetic energy (µJ)')
ax.legend(loc='center right')
ax.text(-.02, 1.02, 'a', transform=ax.transAxes, va='bottom', fontweight='bold', fontsize=9, clip_on=False)

ax = axes[1]
labels = ['Translation\ngain', 'Rotation\nchange', 'Applied\nwork', 'Contact\nnet work']
values = np.array([
    event.Delta_Ktrans_J, event.Delta_Krot_J, event.Wapplied_J, event.Wcontact_net_J]) * 1e6
colors = ['#0072B2', '#6B7280', '#009E73', '#D55E00']
bars = ax.bar(np.arange(4), values, color=colors, width=.68)
ax.axhline(0, color='#111827', lw=.8)
ax.set_xticks(np.arange(4), labels)
ax.set_ylabel('Energy change over impact (µJ)')
pad = .06 * (values.max() - values.min())
for bar, value in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width() / 2, value + (pad if value >= 0 else -pad),
            f'{value:+.3f}', ha='center', va='bottom' if value >= 0 else 'top', fontsize=7)
ax.set_ylim(values.min() - 3 * pad, values.max() + 3 * pad)
ax.text(-.02, 1.02, 'b', transform=ax.transAxes, va='bottom', fontweight='bold', fontsize=9, clip_on=False)

fig.suptitle('First wall impact redirects rotational energy into translation while dissipating net energy', fontsize=9, y=.995)
fig.tight_layout(rect=(0, 0, 1, .96), w_pad=2.0)
fig.savefig(HERE / 'first_impact_energy_transfer.svg', bbox_inches='tight')
fig.savefig(HERE / 'first_impact_energy_transfer.pdf', bbox_inches='tight')
fig.savefig(HERE / 'first_impact_energy_transfer.png', dpi=300, bbox_inches='tight')
plt.close(fig)
