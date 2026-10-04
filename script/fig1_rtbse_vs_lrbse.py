"""
Main figure 1 -- RT-BSE reproduces LR-BSE.

LR-BSE diagonalises the Bethe-Salpeter Hamiltonian; RT-BSE propagates the
density matrix after a delta kick and Fourier transforms the dipole.  They are
two routes to the same spectrum.  Shown here for L = 16 (296 atoms), the
longest ribbon where diagonalisation is still affordable, in the basis the
whole length sweep uses.  Both are TDA with the same broadening, eta = 0.05 eV: a Lorentzian
put in by hand in LR-BSE, the damping applied before the Fourier transform in
RT-BSE (DAMPING 13.1642 fs = hbar / 0.05 eV).

The inset is the ribbon itself, drawn straight from the xyz file.  The kick,
and the polarizability component plotted, are along its long axis.

Figure S2 shows the same comparison for L = 2, 4, 8 and 16 together.

Data: data/16/{lrbse/szv2,rtbse/szv2,struc.xyz}
"""

from plot_param import *
from matplotlib.collections import PolyCollection, PatchCollection
from matplotlib.patches import Circle

L = 16
EMIN, EMAX = 1.0, 8.0

C_LR, C_RT = '#d62728', C_E1

# Ball-and-stick colours and radii (Angstrom) for the structure inset.  The
# stick comes out about a quarter as wide as a carbon ball, which is roughly
# what molecular viewers use by default.
ATOM = {"C": ('#4d4d4d', 0.40), "H": ('#f4f4f4', 0.24)}
COVALENT = {"C": 0.76, "H": 0.31}     # covalent radii [Angstrom]
BOND_TOL = 1.25          # two atoms are bonded within this times r_a + r_b
R_STICK = 0.105          # cylinder radius [Angstrom]
N_STRIP = 11             # strips each cylinder is shaded with
# A fairly oblique light, so the highlight sits up and to the left rather than
# in the middle of each ball, and a soft specular exponent.  A tight highlight
# (exponent ~30) collapses to one or two bright pixels once an atom is only
# ~20 px across, which reads as speckle rather than gloss.
LIGHT = (-0.60, 0.65, 0.45)
SPEC_EXP, SPEC_AMP = 10.0, 0.40
N_SHELL = 32             # nested circles used to shade each sphere
# Where on each nested circle the shading is sampled, as a fraction of its
# radius towards the dark side.  Sampling the darkest point (1.0) under-lights
# the ball by about 40 %; 0.40 reproduces the mean brightness of a per-pixel
# shaded sphere to 3 % for both the carbon and the hydrogen colour.
SAMPLE_W = 0.40


def read_xyz(path):
    """Return (symbols, coordinates) from a plain xyz file."""
    sym, xyz = [], []
    for line in open(path).readlines()[2:]:
        p = line.split()
        if len(p) == 4:
            sym.append(p[0])
            xyz.append([float(v) for v in p[1:]])
    return np.array(sym), np.array(xyz)


def shade(base, nx, ny, nz):
    """Ambient + Lambertian + specular colour for a surface normal (nx, ny, nz).

    One light, placed up and to the left; the viewer looks down +z.  Used for
    both the spheres and the bond cylinders, so the two are lit consistently.
    """
    lgt = np.array(LIGHT) / np.linalg.norm(LIGHT)
    half = lgt + np.array([0.0, 0.0, 1.0])            # light + viewer direction
    half /= np.linalg.norm(half)

    diffuse = np.clip(nx * lgt[0] + ny * lgt[1] + nz * lgt[2], 0.0, 1.0)
    specular = np.clip(nx * half[0] + ny * half[1] + nz * half[2],
                       0.0, 1.0) ** SPEC_EXP
    rgb = np.array(mpl.colors.to_rgb(base))
    out = (rgb * (0.26 + 0.74 * diffuse)[..., None]
           + SPEC_AMP * specular[..., None])
    return np.clip(out, 0.0, 1.0)


def sphere_patches(x, y, radius, colour):
    """A lit sphere drawn as nested circles; returns (patches, facecolours).

    Rasterising each atom and letting matplotlib scale the image down leaves a
    soft, muddy rim, because the sprite is resampled to whatever size the atom
    ends up on the page.  Nested filled circles are vector instead, so the
    outline stays crisp at any size and in the PDF at any zoom.

    The circles shrink from the full radius to nothing while their centres
    drift towards the highlight, which is where a radial gradient gets its
    off-centre look from.  Each is coloured by the true shading sampled a
    fraction SAMPLE_W of its radius towards the dark side, so the outermost
    carries a rim colour and the innermost the specular highlight.
    """
    lgt = np.array(LIGHT) / np.linalg.norm(LIGHT)
    half = lgt + np.array([0.0, 0.0, 1.0])
    half /= np.linalg.norm(half)
    h2 = half[:2]                      # where the highlight sits on the disc
    f = np.linalg.norm(h2)
    d = h2 / f                         # direction from centre to highlight

    patches, colours = [], []
    for k in range(N_SHELL):
        t = k / (N_SHELL - 1.0)
        r = radius * (1.0 - t)
        c = np.array([x, y]) + d * f * radius * t
        u = (c - np.array([x, y])) / radius - d * (r / radius) * SAMPLE_W
        nrm = np.linalg.norm(u)
        if nrm > 1.0:
            u = u / nrm
            nrm = 1.0
        nz = np.sqrt(max(0.0, 1.0 - nrm ** 2))
        patches.append(Circle(c, r))
        colours.append(shade(colour, u[0], u[1], nz))
    return patches, colours


def cylinder_strips(start, end, colour):
    """Shaded cylinder from `start` to `end`, as a list of (quad, rgb) strips.

    A cylinder seen side-on is shaded only across its width: a fraction u of
    the way from one edge to the other, the surface normal is
    u * (perpendicular) + sqrt(1 - u^2) * (towards the viewer).  Splitting the
    bond into strips across its width and colouring each by that normal gives
    the usual rounded look.  Strips overlap a little so that antialiasing
    leaves no seams between them.
    """
    d = end - start
    d = d / np.linalg.norm(d)
    perp = np.array([-d[1], d[0]])

    out = []
    for k in range(N_STRIP):
        u0, u1 = -1 + 2 * k / N_STRIP, -1 + 2 * (k + 1) / N_STRIP
        pad = 0.2 * (u1 - u0)
        a, b = (u0 - pad) * R_STICK, (u1 + pad) * R_STICK
        quad = np.array([start + a * perp, end + a * perp,
                         end + b * perp, start + b * perp])
        u = 0.5 * (u0 + u1)
        out.append((quad, shade(colour, u * perp[0], u * perp[1],
                                np.sqrt(1 - u ** 2))))
    return out


def draw_molecule(ax, sym, xyz):
    """Ball-and-stick picture of a planar molecule, projected on its own plane.

    The ribbon lies in x-y with a z spread of only the sp3 CH2 caps, so the two
    widest axes are the ones worth drawing: the long axis goes horizontal, and
    the narrow one is used as depth.  Atoms are stamped on after the bonds, far
    ones first, so the nearer ones overlap them.

    Two atoms count as bonded within BOND_TOL times the sum of their covalent
    radii.  A single flat distance cutoff does not work here: the two hydrogens
    of each sp3 CH2 cap are only 1.55 A apart, so anything loose enough to
    catch C-C at 1.4 A also draws a bond between them.
    """
    wide = np.argsort(-np.ptp(xyz, axis=0))           # widest axis first
    u, v, depth = (xyz[:, wide[0]], xyz[:, wide[1]], xyz[:, wide[2]])
    pos = np.column_stack([u, v])

    # bonds: each half takes the colour of the atom it grows out of, the way
    # most molecular viewers draw them
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    r = np.array([COVALENT[s] for s in sym])
    bonded = d < BOND_TOL * (r[:, None] + r[None, :])

    quads, rgbs = [], []
    for a, b in zip(*np.where(np.triu(bonded, k=1))):
        mid = 0.5 * (pos[a] + pos[b])
        for k in (a, b):
            for quad, rgb in cylinder_strips(mid, pos[k], ATOM[sym[k]][0]):
                quads.append(quad)
                rgbs.append(rgb)
    ax.add_collection(PolyCollection(quads, facecolors=rgbs, linewidths=0,
                                     zorder=1))

    # atoms, far ones first; one collection keeps the draw order of its patches
    circles, colours = [], []
    for k in np.argsort(depth):                       # back to front
        col, rad = ATOM[sym[k]]
        pk, ck = sphere_patches(u[k], v[k], rad, col)
        circles += pk
        colours += ck
    ax.add_collection(PatchCollection(circles, facecolors=colours,
                                      linewidths=0, match_original=False,
                                      zorder=2))

    ax.set_xlim(u.min() - 0.6, u.max() + 0.6)
    ax.set_ylim(v.min() - 0.6, v.max() + 0.6)
    ax.set_aspect('equal')
    ax.axis('off')


E_lr, y_lr = read_lr(f"{data_dir}/{L}/lrbse/szv2/BSE-TDA-eta=0.050.spectrum")
E_rt, y_rt = read_rt(f"{data_dir}/{L}/rtbse/szv2/"
                     "rec_4-POLARIZABILITY-1_PADE_SPIN_A.dat")
m_lr = (E_lr >= EMIN) & (E_lr <= EMAX)
m_rt = (E_rt >= EMIN) & (E_rt <= EMAX)

# sized for a single REVTeX column: at \columnwidth the tick labels
# come out near 9 pt
fig, ax = plt.subplots(figsize=(9, 6))

ax.fill_between(E_lr[m_lr], 0, y_lr[m_lr], color=tint(C_LR, 0.85), zorder=1)
ax.plot(E_lr[m_lr], y_lr[m_lr], '-', color=C_LR, lw=3.2, zorder=2,
        label="LR-BSE")
ax.plot(E_rt[m_rt], y_rt[m_rt], '--', color=C_RT, lw=2.4, zorder=3,
        label="RT-BSE")

ax.set_xlim(EMIN, EMAX)
ax.set_ylim(0, 1.12 * max(y_lr[m_lr].max(), y_rt[m_rt].max()))
ax.set_xlabel("Energy (eV)")
ax.set_ylabel(r"Im $\alpha_{yy}$  (a.u.)")
ax.grid(True, ls='--', alpha=0.3)
ax.legend(loc='upper right', fontsize=22)

# the ribbon itself, in the empty space under the legend
ins = ax.inset_axes([0.29, 0.49, 0.69, 0.18])
draw_molecule(ins, *read_xyz(f"{data_dir}/{L}/struc.xyz"))

save(fig, "fig1_rtbse_vs_lrbse")
