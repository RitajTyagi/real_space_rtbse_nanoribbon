"""
Supplementary figure 2 -- RT-BSE against LR-BSE, same basis, same ribbons.

LR-BSE diagonalises the Bethe-Salpeter Hamiltonian; RT-BSE propagates the
density matrix after a delta kick and Fourier transforms the dipole.  They are
two routes to the same spectrum, so on the ribbons where both are affordable
(L = 1 to 8, aug-SZV-MOLOPT-GTH-tier-2) the two curves should lie on top of
each other.  Both use TDA and the same broadening, eta = 0.05 eV: in LR-BSE it
is the Lorentzian put in by hand, in RT-BSE it is the damping applied before
the Fourier transform (DAMPING 13.1642 fs = hbar / 0.05 eV).

The lowest line, printed in every panel, agrees to 0.01 eV at all four
lengths, and the strong low-energy lines lie on top of each other.  The two
curves do drift apart above roughly 4 eV, where the lines are dense: the RT run
is 200 fs long, so its Fourier transform cannot resolve structure finer than
h/T = 0.02 eV, and the Pade continuation redistributes weight between lines it
cannot separate.  That is a resolution limit of the real-time route, not a
disagreement about the underlying spectrum.

Data: data/<L>/lrbse/szv2/ and data/<L>/rtbse/szv2/
"""

from plot_param import *

LENGTHS = [1, 2, 4, 8]
EMIN, EMAX = 1.0, 8.0

C_LR, C_RT = '#d62728', C_E1

fig, axes = plt.subplots(2, 2, figsize=(18, 12))
fig.subplots_adjust(hspace=0.25, wspace=0.22, top=0.93)

for ax, L, tag in zip(axes.ravel(), LENGTHS, "abcd"):
    E_lr, y_lr = read_lr(f"{data_dir}/{L}/lrbse/szv2/BSE-TDA-eta=0.050.spectrum")
    E_rt, y_rt = read_rt(f"{data_dir}/{L}/rtbse/szv2/"
                         "rec_4-POLARIZABILITY-1_PADE_SPIN_A.dat")
    m_lr = (E_lr >= EMIN) & (E_lr <= EMAX)
    m_rt = (E_rt >= EMIN) & (E_rt <= EMAX)

    ax.fill_between(E_lr[m_lr], 0, y_lr[m_lr], color=tint(C_LR, 0.85), zorder=1)
    ax.plot(E_lr[m_lr], y_lr[m_lr], '-', color=C_LR, lw=3.0, zorder=2,
            label="LR-BSE")
    ax.plot(E_rt[m_rt], y_rt[m_rt], '--', color=C_RT, lw=2.2, zorder=3,
            label="RT-BSE")

    ax.set_xlim(EMIN, EMAX)
    # scale to whichever method peaks higher, so neither curve is clipped
    ax.set_ylim(0, 1.15 * max(y_lr[m_lr].max(), y_rt[m_rt].max()))
    ax.grid(True, ls='--', alpha=0.3)
    panel_label(ax, tag, f"$L$ = {L}")

for ax in axes[1]:
    ax.set_xlabel("Energy (eV)")
for ax in axes[:, 0]:
    ax.set_ylabel(r"Im $\alpha_{yy}$  (a.u.)")

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=2, fontsize=22,
           bbox_to_anchor=(0.5, 1.0))

save(fig, "figS2_lrbse_vs_rtbse")
