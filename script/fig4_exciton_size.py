"""
Main figure 4 -- how the exciton binding energy and the exciton size track
each other as the ribbon is lengthened.

Left axis: the binding energy of the lowest exciton, E_b = E_g^GW - E1, from
the RT-BSE sweep, available for every length.

Right axis: the size of that same exciton, the descriptor
d_exc = sqrt( <|r_h - r_e|^2> ) printed by the LR-BSE runs, available for
L = 2 to 16, the lengths diagonalisation was affordable for.  The excitation
level used is the one carrying the E1 peak, i.e. the state with the largest y
transition moment within 0.15 eV of the lowest line of the spectrum (y is the
kick direction).

The two are mirror images: the exciton widens from 4.8 to 8.6 A as the ribbon
grows, and stops widening by L = 16, at which point the binding energy also
stops falling.  Beyond that length the ribbon ends no longer squeeze it.

Data: data/csv/exciton_radius.csv and data/csv/gaps_and_excitons.csv.
"""

from plot_param import *

PLOT_FROM = 2                      # L = 1 is a molecule with no peak below 5 eV
C_BIND, C_SIZE = '#2ca02c', '#9467bd'

gaps = pd.read_csv(f"{data_dir}/csv/gaps_and_excitons.csv")
exc = pd.read_csv(f"{data_dir}/csv/exciton_radius.csv")
gaps = gaps[gaps.L >= PLOT_FROM]

fig, ax = plt.subplots(figsize=(9, 6.5))
ax2 = ax.twinx()

l1, = ax.plot(gaps.inv_L, gaps.E_bind_eV, '^-', color=C_BIND, ms=11, lw=2.5,
              label=r"$E_\mathrm{b}$  (RT-BSE)")
l2, = ax2.plot(exc.inv_L, exc.d_exc_A, 's--', color=C_SIZE, ms=11, lw=2.5,
               label=r"$d_\mathrm{exc}$  (LR-BSE)")

ax.set_xlim(-0.02, 0.54)
ax.set_xlabel(r"$1/L$")
ax.set_ylabel(r"Binding energy  $E_\mathrm{b}$ (eV)", color=C_BIND)
ax2.set_ylabel(r"Exciton size  $d_\mathrm{exc}$ (Å)", color=C_SIZE)

# colour each axis to match its curve, so there is no doubt which is which
ax.tick_params(axis='y', colors=C_BIND)
ax2.tick_params(axis='y', colors=C_SIZE)
ax.spines['left'].set_color(C_BIND)
ax2.spines['left'].set_color(C_BIND)
ax.spines['right'].set_color(C_SIZE)
ax2.spines['right'].set_color(C_SIZE)
# the twin axis draws its own frame on top, so turn off the duplicated ticks
ax.tick_params(axis='y', which='both', right=False)
ax2.tick_params(axis='y', which='both', left=False)

ax.grid(True, ls='--', alpha=0.35)
# the band between the two curves at small 1/L is the only clear space
ax.legend(handles=[l1, l2], loc='center left', fontsize=20)

save(fig, "fig4_exciton_size")
