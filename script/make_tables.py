"""
Extract every number the figures need from the raw CP2K output and write it to
data/csv/.  Run this once before the plotting scripts:

    python3 make_tables.py

It produces four tables:

  gaps_and_excitons.csv  PBE gap, G0W0 gap, E1 and the binding energy per L
  peak_tracks.csv        the three tracked peaks E1, E2, E3 per L
  timing.csv             CPU hours for GW, one RT-BSE step and the whole
                         propagation, per L
  exciton_radius.csv     size of the E1 exciton, from the LR-BSE descriptors
"""

from plot_param import *
import re

LENGTHS = [1, 2, 4, 8, 16, 32, 64, 128, 256]

# The three peaks are followed starting from the converged L = 256 spectrum and
# chained back down the sweep.  They cannot all be followed the whole way:
#   E1, E2  are clear from L = 2 on
#   E3      is the intense line near 3.8 eV; below L = 16 the 3.5-4.0 eV region
#           is a dense forest of similar peaks and the assignment would be a
#           coin toss, so the track simply starts at L = 16
TRACK_START = {"E1": 1.56, "E2": 2.07, "E3": 3.82}
TRACK_MIN_L = {"E1": 2, "E2": 2, "E3": 16}
WINDOW = (1.0, 5.0)        # energy window the peaks are searched in [eV]
WINDOW_LR = (1.0, 8.0)     # wider window for the LR-BSE spectra
MATCH_TOL = 0.6            # how far a peak may sit from its predicted position [eV]


def rt_file(L):
    return f"{data_dir}/{L}/rtbse/szv2/rec_4-POLARIZABILITY-1_PADE_SPIN_A.dat"


def grep_value(path, pattern):
    """Return the last number on the first line of `path` containing `pattern`."""
    for line in open(path, errors="ignore"):
        if pattern in line:
            return float(line.split()[-1])
    return np.nan


########## 1. peak tracks ##########
# Collect the candidate peaks of every length first.
peaks = {L: find_peaks(*read_rt(rt_file(L)), *WINDOW) for L in LENGTHS}

# Walk from the longest ribbon down.  Each peak moves up in energy as the
# ribbon gets shorter, and it does so roughly like E(L) = E_inf + A/L, so the
# step it takes grows by a factor of two on every halving of L.  Matching a
# track to the *nearest* peak therefore fails at the short end: between L = 4
# and L = 2 the true step is over 1 eV and the nearest peak is the wrong one.
# Instead each track is extended by extrapolating its last two points linearly
# in 1/L and then taking the peak closest to that prediction.  A track may not
# take a peak another track has already claimed, which is what keeps E1 and E2
# apart at L = 2 (E1 lands on the energy E2 had at L = 4).
tracks = {k: {} for k in TRACK_START}
history = {k: [] for k in TRACK_START}           # list of (1/L, E), longest first
for L in sorted(LENGTHS, reverse=True):
    used = []
    for k in ("E1", "E2", "E3"):
        if L < TRACK_MIN_L[k]:
            continue
        h = history[k]
        if len(h) == 0:
            predict = TRACK_START[k]             # anchor: the L = 256 spectrum
        elif len(h) == 1:
            predict = h[-1][1]
        else:
            (x0, e0), (x1, e1) = h[-2], h[-1]
            predict = e1 + (e1 - e0) / (x1 - x0) * (1.0 / L - x1)
        free = [(E, a) for E, a in peaks[L] if E not in used]
        if not free:
            continue
        E, a = min(free, key=lambda t: abs(t[0] - predict))
        if abs(E - predict) > MATCH_TOL:
            continue
        tracks[k][L] = (E, a)
        history[k].append((1.0 / L, E))
        used.append(E)

rows = []
for L in LENGTHS:
    row = {"L": L, "n_atoms": 18 * L + 8}
    for k in ("E1", "E2", "E3"):
        E, a = tracks[k].get(L, (np.nan, np.nan))
        row[f"{k}_eV"] = E
        row[f"Im_alpha_{k}_per_cell"] = a / L      # per repeat unit
    rows.append(row)
peak_df = pd.DataFrame(rows)
peak_df.to_csv(f"{data_dir}/csv/peak_tracks.csv", index=False, float_format="%.4f")

########## 2. gaps and binding energy ##########
rows = []
for L in LENGTHS:
    out = f"{data_dir}/{L}/rtbse/szv2/output.out"
    pbe = grep_value(out, "SCF indirect band gap (eV):")
    gw = grep_value(out, "G0W0 indirect band gap (eV):")
    e1 = peak_df.loc[peak_df.L == L, "E1_eV"].item()
    rows.append({"L": L, "n_atoms": 18 * L + 8, "inv_L": 1.0 / L,
                 "gap_PBE_eV": pbe, "gap_GW_eV": gw, "E1_eV": e1,
                 "sigma_shift_eV": gw - pbe, "E_bind_eV": gw - e1})
gap_df = pd.DataFrame(rows)
gap_df.to_csv(f"{data_dir}/csv/gaps_and_excitons.csv", index=False,
              float_format="%.4f")

########## 3. timing ##########
# Timings come from the Noctua2 reruns under <L>/rtbse-noctua-timing/.
# Every job used 8 MPI ranks x 16 OpenMP threads per node = 128 cores per node,
# so CPU hours = n_nodes * 128 * seconds / 3600.
CORES_PER_NODE = 8 * 16

rows = []
for L in LENGTHS:
    out = f"{data_dir}/{L}/rtbse-noctua-timing/output.out"
    if not os.path.exists(out):
        continue                                   # L = 1 was not rerun there

    # node count from the batch script
    nodes = int([l.split("-N")[1] for l in open(f"{data_dir}/{L}/rtbse-noctua-timing/submit.slurm")
                 if l.startswith("#SBATCH -N")][0])
    cores = nodes * CORES_PER_NODE

    lines = open(out, errors="ignore").readlines()

    # (a) per-step RT-BSE wall time, from the "RTBSE| <step> <time>" table.
    # The first ten steps still include start-up effects and are dropped; the
    # mean of the rest is the quoted cost.  Taking the median instead moves it
    # by under 1 % everywhere except L = 32, where a tail of slow steps pulls
    # the mean up by 3.6 %, and it does not change the fitted exponents.
    steps = np.array([float(l.split()[2]) for l in lines
                      if re.match(r"\s*RTBSE\|\s+\d+\s+[\d.]+", l)])
    step_s = float(np.mean(steps[10:]))
    # Total propagation cost: the step times actually printed, added up.  No
    # step count is assumed -- each run propagated for as long as it ran.  This
    # sum reproduces CP2K's own solve_rk4_timestep timer to better than 1 %.
    total_s = float(steps.sum())

    # (b) total GW wall time.  Finished runs print a "gw" row in the closing
    # timing report.  L = 256 is still propagating, so its report does not
    # exist yet and the GW cost is rebuilt by adding up the stage timers that
    # the GW section prints as it goes (grid optimisation, Z_lP, chi, W, Sigma).
    # On the finished runs that sum recovers 95-99 % of the reported total, so
    # it is divided by the L = 128 ratio to put it on the same footing.
    gw_row = [l for l in lines if l.split()[:2] == ["gw", "1"]]
    stage = sum(float(l.split("time")[-1].replace(":", "").split()[0])
                for l in lines
                if ("Execution time" in l or "execution time" in l)
                and "for atom" not in l
                and not any(s in l for s in ("DOS, LDOS", "V_aux", "W0_grid")))
    if gw_row:
        gw_s, gw_complete = float(gw_row[0].split()[-1]), True
    else:
        gw_s, gw_complete = stage / 0.950, False

    rows.append({"L": L, "n_atoms": 18 * L + 8, "nodes": nodes, "cores": cores,
                 "n_steps_done": len(steps), "gw_seconds": gw_s,
                 "gw_cpuh": cores * gw_s / 3600.0,
                 "rtbse_step_seconds": step_s,
                 "rtbse_step_cpuh": cores * step_s / 3600.0,
                 "rtbse_total_seconds": total_s,
                 "rtbse_total_cpuh": cores * total_s / 3600.0,
                 "gw_complete": gw_complete,
                 "stage_sum_seconds": stage})
time_df = pd.DataFrame(rows)
time_df.to_csv(f"{data_dir}/csv/timing.csv", index=False, float_format="%.4f")

########## 4. size of the E1 exciton ##########
# The LR-BSE runs print exciton descriptors per excitation level.  To pick the
# level that makes the E1 peak, take the one with the largest y transition
# moment within 0.15 eV of the lowest line of the spectrum -- y is the kick
# direction, so d_y^2 is what the plotted spectrum weighs each state by.
LR_LENGTHS = [2, 4, 8, 16]       # the lengths LR-BSE was affordable for
MATCH_WINDOW = 0.15              # eV


def parse_lrbse(path):
    """Excitation energies, y transition moments and exciton descriptors.

    CP2K prints the isotropic descriptor table and then a per-direction one
    whose rows have the same shape, so parsing has to stop at the second
    header or the first table is silently overwritten.
    """
    E, dy, dexc, dexc_y = {}, {}, {}, {}
    sec = None
    for line in open(path, errors="ignore"):
        if "Excitation energies from solving the BSE" in line: sec = "E"; continue
        if "Optical properties from solving the BSE" in line: sec = "D"; continue
        if "Exciton descriptors per direction" in line: sec = "XD"; continue
        if "Exciton descriptors from solving the BSE" in line: sec = "X"; continue
        if not line.startswith(" BSE|"):
            continue
        p = line.split()[1:]
        if not (p and p[0].isdigit()):
            continue
        if sec == "E" and len(p) == 4:
            E[int(p[0])] = float(p[3])
        elif sec == "D" and len(p) == 6:
            dy[int(p[0])] = float(p[3])
        elif sec == "X" and len(p) == 7:
            dexc[int(p[0])] = float(p[5])
        elif sec == "XD" and len(p) >= 6 and p[1] == "y":
            dexc_y[int(p[0])] = float(p[5])
    return E, dy, dexc, dexc_y


rows = []
for L in LR_LENGTHS:
    base = f"{data_dir}/{L}/lrbse/szv2"
    E, dy, dexc, dexc_y = parse_lrbse(f"{base}/output.out")
    Es, ys = read_lr(f"{base}/BSE-TDA-eta=0.050.spectrum")
    peak = min(find_peaks(Es, ys, *WINDOW_LR))[0]
    near = {n: dy[n] ** 2 for n in E
            if abs(E[n] - peak) < MATCH_WINDOW and n in dy}
    n1 = max(near, key=near.get)
    rows.append({"L": L, "inv_L": 1.0 / L, "n_exc": n1,
                 "omega_eV": E[n1], "dy2": dy[n1] ** 2,
                 "d_exc_A": dexc[n1], "d_exc_y_A": dexc_y[n1],
                 "E_bind_eV": gap_df.loc[gap_df.L == L, "E_bind_eV"].item()})
exc_df = pd.DataFrame(rows)
exc_df.to_csv(f"{data_dir}/csv/exciton_radius.csv", index=False,
              float_format="%.4f")

print(peak_df.to_string(index=False), "\n")
print(gap_df.to_string(index=False), "\n")
print(time_df.to_string(index=False), "\n")
print(exc_df.to_string(index=False))
