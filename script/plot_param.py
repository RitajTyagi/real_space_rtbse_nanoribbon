"""
Shared matplotlib settings and small file readers for every figure in this repo.

Every figure script starts with  `from plot_param import *`  so that all panels
share one look (serif/STIX fonts, inward ticks, same colours).

Units note: CP2K writes the RT-BSE polarizability on a frequency grid in
Hartree, the LR-BSE spectrum on a grid already in eV.  Both files store the
same quantity -- Im alpha_mu,mu'(omega) in atomic units -- so once the RT axis
is converted the two can be plotted on top of each other without rescaling.
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import os

########## Paths ##########
# work_dir is the folder holding this file, i.e. <repo>/script
work_dir = sys.path[0]
data_dir = os.path.join(work_dir, "..", "data")
fig_dir = os.path.join(work_dir, "..", "figure")

########## Physical constants ##########
HARTREE_EV = 27.211386245988          # 1 Hartree in eV

########## Plotting parameters ##########
plt.rcParams['axes.labelsize'] = 26
plt.rcParams['axes.titlesize'] = 24
plt.rcParams.update({'font.size': 22})
plt.rcParams['axes.linewidth'] = 1.6
plt.rcParams['xtick.major.pad'] = 6
plt.rcParams['ytick.major.pad'] = 6
plt.rcParams['xtick.major.size'] = 7
plt.rcParams['ytick.major.size'] = 7
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['xtick.major.width'] = 1.6
plt.rcParams['ytick.major.width'] = 1.6
plt.rcParams['axes.labelpad'] = 12
plt.rcParams['lines.linewidth'] = 2.5
plt.rcParams['legend.frameon'] = False
plt.rcParams['figure.figsize'] = [15, 10]
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.serif"] = ["STIXGeneral", "Times New Roman", "DejaVu Serif"]
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
# keep text as text (not outlines) in the PDF, so it stays editable/searchable
mpl.rcParams['pdf.fonttype'] = 42

########## Colours ##########
colors = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']

C_E1, C_E2, C_E3 = '#1f77b4', '#ff7f0e', '#2ca02c'   # the three tracked peaks
C_GREY, C_GREY_LINE = '#d5d4cd', '#8a8a84'           # ghosted reference curve
INK_GREY = '#55544f'                                 # secondary annotation text
INK = '#111111'                                      # panel labels, primary text

# Blue ramp used for the stacked spectra, light (short ribbon) to dark (long).
# Five anchors is the longest single-hue ramp that still reads as distinct
# steps; intermediate rows are interpolated between them.
BLUE_RAMP = ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#104281']

########## File readers ##########

def read_rt(path):
    """Read a CP2K RT-BSE *POLARIZABILITY*.dat file.

    Columns are: omega [Hartree], then (Re, Im) pairs for alpha_xy, alpha_yy,
    alpha_zy -- the second index is the kick direction, which is y (the ribbon
    axis) in every run here.  We return the yy component, i.e. the response
    along the kick.
    """
    d = np.loadtxt(path)
    return d[:, 0] * HARTREE_EV, d[:, 4]


def read_lr(path):
    """Read a CP2K LR-BSE *.spectrum file and return (E [eV], Im alpha_yy).

    Columns are: frequency [eV], the average, then the nine Cartesian
    components in the order xx xy xz yx yy yz zx zy zz -- so yy is column 6.
    """
    d = np.loadtxt(path)
    return d[:, 0], d[:, 6]


def find_peaks(E, y, emin, emax, rel_floor=0.02):
    """Return the local maxima of y(E) inside [emin, emax] as a list of (E, y).

    Only maxima taller than `rel_floor` times the tallest peak in the window
    are kept, which removes the small numerical ripples of the Pade
    continuation without touching any real line.
    """
    m = (E >= emin) & (E <= emax)
    E, y = E[m], y[m]
    i = np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    if len(i) == 0:
        return []
    floor = rel_floor * y[i].max()
    return [(E[k], y[k]) for k in i if y[k] >= floor]


def tint(color, frac):
    """Blend `color` `frac` of the way towards white and return it opaque.

    The figures are saved with a transparent background, and a semi-transparent
    patch drawn on a transparent background keeps its full saturation in the
    file -- it then looks pale on a white page but garish on a dark slide.
    Pre-blending with white here instead of passing alpha= keeps every fill
    looking the same whatever it is placed on.
    """
    rgb = np.array(mpl.colors.to_rgb(color))
    return tuple(rgb + frac * (1.0 - rgb))


def ramp_color(frac):
    """Colour at position `frac` in [0, 1] along BLUE_RAMP (linear in sRGB)."""
    anchors = np.array([mpl.colors.to_rgb(c) for c in BLUE_RAMP])
    x = frac * (len(anchors) - 1)
    lo = int(np.clip(np.floor(x), 0, len(anchors) - 2))
    t = x - lo
    return tuple((1 - t) * anchors[lo] + t * anchors[lo + 1])


def panel_label(ax, letter, extra="", fontsize=24, newline=False,
                y=0.955):
    """Put "(a)" in the top-left corner INSIDE the axes.

    Panels carry only their letter; what they show belongs in the caption.
    `extra` adds a short identifying label where the panel would otherwise be
    ambiguous (the ribbon length, a keyword name).  `newline` puts that label
    on its own line, which keeps a long keyword from running into a legend.
    """
    sep = "\n" if newline else "  "
    ax.text(0.025, y, f"({letter})" + (sep + extra if extra else ""),
            transform=ax.transAxes, ha='left', va='top',
            fontsize=fontsize, color=INK, linespacing=1.4)


def save(fig, name):
    """Save a figure as transparent PDF and PNG into <repo>/figure."""
    for ext, kw in (("pdf", {}), ("png", {"dpi": 400})):
        fig.savefig(os.path.join(fig_dir, f"{name}.{ext}"),
                    bbox_inches='tight', pad_inches=0.05,
                    transparent=True, **kw)
    print(f"wrote figure/{name}.pdf and figure/{name}.png")
