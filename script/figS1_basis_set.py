"""
Supplementary figure 1 -- basis-set check on the LR-BSE spectra.

The length sweep uses aug-SZV-MOLOPT-GTH-tier-2, a small basis.  Here the same
LR-BSE spectrum is recomputed with aug-cc-pVDZ for the two ribbons where that
is still affordable and the low-energy region is already formed.

Over 0-5 eV the two bases agree on every strong line to about 0.1 eV in
position, with heights differing by up to ~20 %.  The lowest line, which is the
one the length sweep follows, is at the same energy in both bases at L = 4.

Data: data/<L>/lrbse/{avdz,szv2}/BSE-TDA-eta=0.050.spectrum
"""

from plot_param import *

LENGTHS = [2, 4]
EMIN, EMAX = 0.0, 5.0

C_AVDZ = '#d62728'
C_SZV = C_E1

fig, axes = plt.subplots(1, 2, figsize=(16, 7))
fig.subplots_adjust(wspace=0.18, top=0.9)

for ax, L, tag in zip(axes, LENGTHS, "ab"):
    E_s, y_s = read_lr(f"{data_dir}/{L}/lrbse/szv2/BSE-TDA-eta=0.050.spectrum")
    E_a, y_a = read_lr(f"{data_dir}/{L}/lrbse/avdz/BSE-TDA-eta=0.050.spectrum")
    ms = (E_s >= EMIN) & (E_s <= EMAX)
    ma = (E_a >= EMIN) & (E_a <= EMAX)

    # the basis the sweep uses is filled, the reference basis goes on top as a
    # line, so neither curve is hidden by the other's fill
    ax.fill_between(E_s[ms], 0, y_s[ms], color=tint(C_SZV, 0.82), zorder=1)
    ax.plot(E_s[ms], y_s[ms], '-', color=C_SZV, lw=2.5, zorder=2,
            label="aug-SZV-MOLOPT-GTH-tier-2")
    ax.plot(E_a[ma], y_a[ma], '--', color=C_AVDZ, lw=2.2, zorder=3,
            label="aug-cc-pVDZ")

    ax.set_xlim(EMIN, EMAX)
    # scale to whichever basis peaks higher, so neither curve is clipped
    ax.set_ylim(0, 1.15 * max(y_s[ms].max(), y_a[ma].max()))
    ax.set_xlabel("Energy (eV)")
    ax.grid(True, ls='--', alpha=0.3)
    panel_label(ax, tag, f"$L$ = {L}")

axes[0].set_ylabel(r"Im $\alpha_{yy}$  (a.u.)")
# one legend for all three panels, above the figure
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=2, fontsize=21,
           bbox_to_anchor=(0.5, 1.0))

save(fig, "figS1_basis_set")
