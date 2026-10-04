"""
Main figure 3 -- cost of the two stages, in CPU hours.

  * the total G0W0 time, once per calculation
  * the total RT-BSE propagation time of the run
  * the mean wall time of a single RT-BSE propagation step

The total is the step times printed by the run, added up -- no step count is
assumed anywhere, each job simply propagated for as long as it ran.  That sum
reproduces CP2K's own solve_rk4_timestep timer to better than 1 %.

Both are converted from wall time to CPU hours with the resources the job
actually asked for: every job ran 8 MPI ranks x 16 OpenMP threads per node,
i.e. 128 cores per node, so CPU hours = n_nodes * 128 * seconds / 3600.

Both axes are logarithmic, so a power law t = A L^p is a straight line and p is
its slope.  The dashed lines are fitted to the three largest ribbons only
(L = 64, 128, 256), which is the asymptotic regime, and run a little past the
last point; the exponent is given in the legend.  Shorter ribbons fall below
those lines because the fixed cost of the RI-RS grid optimisation still
dominates there.  The ribbon length L is used rather than the atom count
because N = 18 L + 8 is proportional to L for all but the shortest ribbons, so
the two give the same exponent to within 0.01.

These timings come from the Noctua2 reruns (data/<L>/rtbse-noctua-timing/),
not from the production runs the spectra come from -- see the README.

Data: data/csv/timing.csv (built by make_tables.py).
"""

from plot_param import *

N_FIT = 3                # number of largest ribbons the power law is fitted to

df = pd.read_csv(f"{data_dir}/csv/timing.csv")

# sized for a single REVTeX column: at \columnwidth the tick labels
# come out near 9 pt
fig, ax = plt.subplots(figsize=(8.5, 7.8))

SERIES = [("rtbse_total_cpuh", "Total RT-BSE", '#2ca02c', '^'),
          ("gw_cpuh", r"Total $G_0W_0$", '#d62728', 's'),
          ("rtbse_step_cpuh", "One RT-BSE step", C_E1, 'o')]

for col, label, c, marker in SERIES:
    L, t = df.L.values, df[col].values

    # The measured points are joined by a thin solid line.  The power law
    # fitted to the last N_FIT points is drawn as a thick dash on top -- the
    # fit is so close to the data there that a thin dashed line would simply
    # disappear under the solid one -- and runs a little past the last point,
    # so that it reads as a trend rather than stopping dead on the marker.
    p, logA = np.polyfit(np.log(L[-N_FIT:]), np.log(t[-N_FIT:]), 1)
    Lf = np.logspace(np.log10(L[-N_FIT]), np.log10(L[-1] * 1.35), 50)

    ax.plot(L, t, '-', color=c, lw=1.8, zorder=1,
            label=f"{label}   " + rf"$\propto L^{{{p:.1f}}}$")
    ax.plot(Lf, np.exp(logA) * Lf ** p, '--', color=c, lw=4.0,
            dashes=(3.5, 2.2), zorder=2)
    ax.plot(L, t, marker, color=c, ms=13, zorder=3)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks(df.L.values)
ax.set_xticklabels([str(L) for L in df.L.values])
ax.minorticks_off()
ax.set_xlabel(r"Repeat units  $L$")
ax.set_ylabel("Computational cost (CPU hours)")
ax.grid(True, which='major', ls='--', alpha=0.35)
ax.legend(loc='upper left', fontsize=19)

# second x axis giving the atom count that goes with each length, N = 18 L + 8
top = ax.secondary_xaxis('top')
top.set_xscale('log')
top.set_xticks(df.L.values)
top.set_xticklabels([str(n) for n in df.n_atoms.values], fontsize=18)
top.minorticks_off()
top.set_xlabel("Number of atoms", labelpad=14)

save(fig, "fig3_timing")
