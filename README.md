# linRT-BSE length sweep of a 7-AGNR nanoribbon

CP2K RI-RS `G0W0` + linear-response real-time BSE on an armchair graphene
nanoribbon of `L` repeat units, `L = 1 ... 256` (26 to 4616 atoms), together
with the LR-BSE and parameter-convergence checks behind it.

## Layout

```
data/                       raw CP2K input and output, one folder per L
  <L>/struc.xyz             geometry
  <L>/rtbse/szv2/           RT-BSE production run, aug-SZV-MOLOPT-GTH-tier-2
  <L>/lrbse/szv2/           LR-BSE, same basis            (L = 1, 2, 4, 8)
  <L>/lrbse/avdz/           LR-BSE, aug-cc-pVDZ           (L = 1, 2, 4)
  <L>/rtbse-noctua-timing/  timing rerun, incl. submit.slurm  (L = 2 ... 256)
  16/conv_parameter/<par>/<value>/   RI-RS / RT-BSE parameter scans
  csv/                      tables extracted by script/make_tables.py

script/                     one script per figure, plus
  plot_param.py             shared matplotlib settings and file readers
  make_tables.py            writes data/csv/*.csv

figure/                     figures 
```

Geometry: `n_atoms = 18 L + 8`, repeat unit 4.264 Å.
Every run uses PBE/GPW, `CUTOFF 600`, `G0W0` with 30 time/frequency points,
TDA, and a 200 fs propagation with `DAMPING 13.1642 fs`, i.e. the same
broadening `eta = 0.05 eV` in RT-BSE and LR-BSE.
Spectra are the Padé-continued `*_PADE_*` files, `Im alpha_yy`, the component
along the kick (the ribbon axis).

## Figures

| file | content | data |
|---|---|---|
| `fig1_rtbse_vs_lrbse` | RT-BSE against LR-BSE at `L = 8`, with the ribbon | `8/{lrbse,rtbse}/szv2`, `8/struc.xyz` |
| `fig2_length_sweep` | (a) stacked spectra per `L`; (b) the three main peaks vs `1/L`; (c) PBE gap, `G0W0` gap, `E1` and binding energy vs `1/L` | `<L>/rtbse/szv2` |
| `fig3_timing` | CPU hours: total GW, and one RT-BSE step | `<L>/rtbse-noctua-timing` |
| `figS1_basis_set` | LR-BSE, aug-cc-pVDZ vs aug-SZV-MOLOPT-GTH-tier-2 | `<L>/lrbse` |
| `figS2_lrbse_vs_rtbse` | RT-BSE against LR-BSE, `L = 1 ... 8` | `<L>/{lrbse,rtbse}/szv2` |
| `figS3_rtbse_parameters` | spectrum vs the four RI-RS / RT-BSE cutoffs | `16/conv_parameter` |


## Note on the timings

The timings in `fig3_timing` do **not** come from the production runs the
spectra come from.  Those ran on a machine where the per-step cost was not
reproducible between jobs, so wall times from them are not comparable across
`L`.  The whole sweep was therefore rerun on Noctua2, where the step time is
stable, and only those reruns are used for timing.  Everything else in this
repository comes from the production runs.

The quoted step cost is the mean over the completed steps, after dropping the
first ten, which still carry start-up effects.  The `L = 256` job was still
propagating when the figures were made, so for it that mean runs over the steps
finished so far.

The dashed power laws in `fig3` are fitted to the last three largest systems only
(`L` = 64, 128, 256).

## Rebuilding

```
cd script
python3 make_tables.py        # writes data/csv/*.csv
python3 fig1_rtbse_vs_lrbse.py  # and the other fig*.py / figS*.py
```

Needs numpy, pandas and matplotlib only.
