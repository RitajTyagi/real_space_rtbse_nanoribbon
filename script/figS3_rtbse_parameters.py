"""
Supplementary figure 3 -- convergence of the L = 16 spectrum with the four
RI-RS / RT-BSE parameters that control the approximations.

  CUTOFF_ATOMIC_CLUSTER   size of the atomic cluster the RI-RS grid is
                          optimised on                               [Angstrom]
  CUTOFF_RADIUS_RL_AO     range kept in the real-space AO product basis
                                                                     [Angstrom]
  CUTOFF_RADIUS_W0        truncation of W in the BSE kernel           [Angstrom]
  N_POINT_PERCENTAGE      fraction of candidate points the grid optimiser keeps
                                                                             [%]

Only two of the four matter here: the AO range has to reach 6 Angstrom and the
grid optimiser has to keep at least 15 % of its points.  Below those the lowest
line moves by 0.21 and 0.02 eV respectively.  The W0 truncation does nothing
between 30 and 60 Angstrom, which matters for the length sweep: the production
runs use 50 Angstrom, so the binding energy plateau in figure 2(d) is not an
artefact of that cutoff.

Each value gets its own colour *and* its own dash pattern, drawn thin, so that
curves which nearly coincide can still be told apart: where they agree the dash
patterns interleave, where they differ one line separates from the bundle.

Data: data/16/conv_parameter/<parameter>/<value>/
"""

from plot_param import *

EMIN, EMAX = 1.0, 5.0

# One colour and one dash pattern per value, so overlapping curves interleave
# instead of hiding each other.  The dash patterns are taken from the END of
# this list, so the largest (most converged) value is always the solid line and
# the looser ones are progressively more broken up.
CURVE_C = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4']
CURVE_LS = [(0, (1, 2)), (0, (7, 2, 1, 2)), (0, (6, 3)), '-']

# directory name, label, values, production value, unit
SETS = [("CUTOFF_ATOMIC_CLUSTER", "CUTOFF_ATOMIC_CLUSTER", [2, 3, 4], 3, "Å"),
        ("CUTOFF_RADIUS_RL_AO", "CUTOFF_RADIUS_RL_AO", [5, 6, 7, 8], 7, "Å"),
        ("CUTOFF_RADIUS_W0", "CUTOFF_RADIUS_W0", [30, 40, 50, 60], 50, "Å"),
        ("GRID_PERCENTAGE", "N_POINT_PERCENTAGE", [10, 15, 20, 25], 20, "%")]

fig, axes = plt.subplots(2, 2, figsize=(18, 12))
fig.subplots_adjust(hspace=0.35, wspace=0.20)

for ax, (folder, label, values, used, unit), tag in zip(axes.ravel(), SETS, "abcd"):
    for i, v in enumerate(values):
        E, y = read_rt(f"{data_dir}/16/conv_parameter/{folder}/{v}/"
                       "rec_4-POLARIZABILITY-1_PADE_SPIN_A.dat")
        m = (E >= EMIN) & (E <= EMAX)
        tail = "  (used)" if v == used else ""
        # zorder decreasing with the value: the solid, most converged curve
        # goes underneath and the broken ones are drawn over it, so where they
        # coincide the dashes still show through instead of being covered
        ax.plot(E[m], y[m], color=CURVE_C[i], lw=1.8, zorder=10 - i,
                ls=CURVE_LS[len(CURVE_LS) - len(values) + i],
                label=f"{v} {unit}{tail}")

    ax.set_xlim(EMIN, EMAX)
    ax.set_ylim(0, 1.22 * y[m].max())        # headroom for the legend
    ax.grid(True, ls='--', alpha=0.3)
    # Keyword and panel letter go in the title above the axes.  Inside the
    # panel the keyword ends up level with the topmost y tick label and the two
    # read as one run of text, and there is no other free space for it.
    ax.set_title(f"({tag})  {label}", loc='left', pad=12, fontsize=22)
    ax.legend(loc='upper right', fontsize=18)

for ax in axes[1]:
    ax.set_xlabel("Energy (eV)")
for ax in axes[:, 0]:
    ax.set_ylabel(r"Im $\alpha_{yy}$  (a.u.)")

save(fig, "figS3_rtbse_parameters")
