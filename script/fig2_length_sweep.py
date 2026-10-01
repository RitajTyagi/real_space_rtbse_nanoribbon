"""
Main figure 2 -- the linRT-BSE length sweep of the 7-AGNR nanoribbon.

(a) stacked absorption spectra, one row per ribbon length, each scaled to its
    own maximum, with the converged L = 256 spectrum ghosted in grey behind
    every row
(b) the three main peaks against 1/L
(c) PBE gap, G0W0 gap, E1 and the binding energy against 1/L

Both (b) and (c) show L >= 16, the range where E(L) = E_inf + A/L holds, so
every quantity is a straight line whose intercept at 1/L = 0 is the value for
the infinite ribbon.  Panel (c) has a broken y axis: the G0W0 gap sits twice as
high as the other three, which would otherwise be squeezed into a tenth of the
panel.

Data: data/<L>/rtbse/szv2/ and the tables in data/csv/.
"""

from plot_param import *

LENGTHS = [2, 4, 8, 16, 32, 64, 128, 256]     # L = 1 has no peak below 5 eV
EMIN, EMAX = 1.0, 5.0                         # energy window shown [eV]
REF = 256                                     # length used as the grey reference
FIT_FROM = 16                                 # 1/L is linear only from here on
XMAX = 0.082                                  # room at the right for the numbers

TRACKS = [("E1", C_E1), ("E2", C_E2), ("E3", C_E3)]

# Where the "E = x.xx eV" label goes in (b) relative to its own line: E1 and E2
# have room above, E3 nearly touches the top of the panel.
LABEL_DY = {"E1": 0.11, "E2": 0.11, "E3": -0.30}

# column, legend label, colour, marker  -- for panel (c)
QUANTITIES = [("gap_GW_eV", r"$G_0W_0$ gap", '#d62728', 's'),
              ("E1_eV", r"$E_1$", C_E1, 'D'),
              ("E_bind_eV", r"$E_\mathrm{b}$", '#2ca02c', '^'),
              ("gap_PBE_eV", "PBE gap", '#9467bd', 'o')]

BAND_HI = (3.11, 3.35)       # (c) upper band: the G0W0 gap
BAND_LO = (1.49, 1.79)       # (c) lower band: E1, E_b and the PBE gap
                             # the top of it is left clear for the legend

peaks = pd.read_csv(f"{data_dir}/csv/peak_tracks.csv")
gaps = pd.read_csv(f"{data_dir}/csv/gaps_and_excitons.csv")


def spectrum(L):
    """Im alpha_yy of ribbon L inside the plotting window, per repeat unit."""
    E, y = read_rt(f"{data_dir}/{L}/rtbse/szv2/"
                   "rec_4-POLARIZABILITY-1_PADE_SPIN_A.dat")
    m = (E >= EMIN) & (E <= EMAX)
    return E[m], y[m] / L


def fit_line(L, v):
    """Least-squares straight line in 1/L; returns (slope, value at 1/L = 0)."""
    slope, v_inf = np.polyfit(1.0 / np.asarray(L), np.asarray(v), 1)
    return slope, v_inf


fig = plt.figure(figsize=(21, 11))
# (c) is given more height than (b) and the gap between them is kept tight:
# (b) is three nearly flat lines that need little room, while (c) has four
# curves packed into a narrow energy range and benefits from every pixel.
gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1], height_ratios=[1, 1.5],
                      hspace=0.17, wspace=0.22)
axL = fig.add_subplot(gs[:, 0])                       # (a)
axB = fig.add_subplot(gs[0, 1])                       # (b)
# (c) is one plot drawn on two stacked axes so that the y axis can be broken
gsC = gs[1, 1].subgridspec(2, 1, height_ratios=[1, 1.9], hspace=0.10)
axC_hi = fig.add_subplot(gsC[0])
axC_lo = fig.add_subplot(gsC[1])

########## (a) stacked spectra ##########
E_ref, y_ref = spectrum(REF)
y_ref = y_ref / y_ref.max()

# light vertical guides at the converged peak positions, so the three lines
# followed in (b) can be picked out in the spectra
for name, c in TRACKS:
    axL.axvline(peaks.loc[peaks.L == REF, f"{name}_eV"].item(),
                color=c, lw=1.2, ls=':', alpha=0.55, zorder=0)

for row, L in enumerate(LENGTHS):
    E, y = spectrum(L)
    base = len(LENGTHS) - 1 - row          # row 0 drawn at the top
    peak_height = y.max()
    y = 0.88 * y / peak_height             # 0.88 leaves a gap between rows
    c = ramp_color(row / (len(LENGTHS) - 1))

    axL.fill_between(E_ref, base, base + 0.88 * y_ref, color=C_GREY,
                     zorder=2 * row + 1)
    axL.fill_between(E, base, base + y, color=tint(c, 0.60), zorder=2 * row + 2)
    # the reference outline goes on top, so the grey shape stays readable where
    # this row's own fill covers it
    axL.plot(E_ref, base + 0.88 * y_ref, color=C_GREY_LINE, lw=1.2,
             zorder=2 * row + 2)
    axL.plot(E, base + y, color=c, lw=2.2, zorder=2 * row + 2)

    axL.text(EMIN + 0.05, base + 0.60, f"L = {L}", color=c,
             fontsize=22, fontweight='bold', va='bottom')

axL.set_xlim(EMIN, EMAX)
axL.set_ylim(-0.05, len(LENGTHS) + 0.65)   # headroom for the panel label
axL.set_yticks([])
axL.set_xlabel("Energy (eV)")
axL.set_ylabel(r"Im $\alpha_{yy}$ per repeat unit")
panel_label(axL, "a", y=0.98)

########## (b) the three main peaks ##########
for name, c in TRACKS:
    d = peaks.dropna(subset=[f"{name}_eV"])
    d = d[d.L >= FIT_FROM]
    L, E = d.L.values, d[f"{name}_eV"].values

    slope, e_inf = fit_line(L, E)
    xf = np.linspace(0, 1.0 / FIT_FROM, 50)
    axB.plot(xf, slope * xf + e_inf, '--', color=c, lw=2.0)
    axB.plot(1.0 / L, E, 'o', color=c, ms=11)
    axB.plot(0, e_inf, 'o', ms=14, mfc='white', mec=c, mew=2.8)
    axB.text(0.0015, e_inf + LABEL_DY[name], f"{name} = {e_inf:.2f} eV",
             color=c, fontsize=20, fontweight='bold')

axB.set_xlim(-0.004, 1.0 / FIT_FROM + 0.004)
axB.set_xticks([0.00, 0.02, 0.04, 0.06])
axB.set_ylim(1.42, 4.30)             # headroom above E3 for the panel label
# pin the ticks: the panel is short, and autoscaling drops to whole eV there
axB.set_yticks([1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
axB.set_ylabel("Peak energy (eV)")
axB.grid(True, ls='--', alpha=0.35)
panel_label(axB, "b")

########## (c) gaps, exciton and binding energy ##########
# every series is drawn on both axes; each axis then clips to its own band
for col, label, c, marker in QUANTITIES:
    d = gaps[gaps.L >= FIT_FROM]
    L, v = d.L.values, d[col].values
    slope, v_inf = fit_line(L, v)
    xf = np.linspace(0, 1.0 / FIT_FROM, 50)

    for ax, (lo, hi) in ((axC_hi, BAND_HI), (axC_lo, BAND_LO)):
        ax.plot(xf, slope * xf + v_inf, '--', color=c, lw=1.8)
        # the label is attached on the lower axis only, where the legend lives;
        # every curve is drawn on both axes and simply clipped by the band
        ax.plot(1.0 / L, v, marker + '-', color=c, ms=10, lw=2.2,
                label=label if ax is axC_lo else None)
        ax.plot(0, v_inf, marker, ms=13, mfc='white', mec=c, mew=2.6)
        # Extrapolated limit, written past the shortest ribbon shown, on the
        # axis whose band actually holds this curve -- text is not clipped by
        # default, so drawing it on both would push it outside the figure and
        # wreck the bounding box.
        if lo <= v[0] <= hi:
            ax.text(1.0 / FIT_FROM + 0.004, v[0], f"{v_inf:.2f}",
                    color=c, fontsize=19, fontweight='bold', va='center')

axC_hi.set_ylim(*BAND_HI)
axC_lo.set_ylim(*BAND_LO)
axC_hi.set_yticks([3.2, 3.3])
axC_lo.set_yticks([1.5, 1.6, 1.7])
for ax in (axC_hi, axC_lo):
    ax.set_xlim(-0.004, XMAX)
    ax.grid(True, ls='--', alpha=0.35)
axC_hi.set_xticklabels([])
axC_lo.set_xticks([0.00, 0.02, 0.04, 0.06])
axC_lo.set_xlabel(r"$1/L$")
axC_lo.set_ylabel("Energy (eV)")
axC_lo.yaxis.set_label_coords(-0.105, 0.75)
# two columns, so the legend is two rows tall and sits inside the empty strip
# left above E1 rather than running down into the curves
axC_lo.legend(loc='upper left', fontsize=16, ncol=2, handlelength=1.6,
              borderpad=0.3, labelspacing=0.3, handletextpad=0.5,
              columnspacing=1.2)

# break the axis: drop the facing spines and draw the usual pair of slashes
axC_hi.spines['bottom'].set_visible(False)
axC_lo.spines['top'].set_visible(False)
# the shared style puts ticks on all four sides; on the two spines that the
# break removes they would be left hanging in the gap, so turn them off
axC_hi.tick_params(bottom=False, labelbottom=False)
axC_lo.tick_params(top=False)
kw = dict(marker=[(-1, -0.6), (1, 0.6)], ms=11, mew=1.6, color=INK,
          ls='none', clip_on=False)
axC_hi.plot([0, 1], [0, 0], transform=axC_hi.transAxes, **kw)
axC_lo.plot([0, 1], [1, 1], transform=axC_lo.transAxes, **kw)
panel_label(axC_hi, "c")

save(fig, "fig2_length_sweep")
