"""
=========================================================================
DOE Analysis - Baseline (B) vs Moderate (M) block
Response: % change in Motility Metric from the healthy baseline

Levels:  B = baseline   eta = 0; vartheta^h = 1.0
         M = moderate   eta = 0.5*eta_ref; vartheta^h = 2.0
                          for the hypertrophy factor

Python port of the original MATLAB script. Console output, figure layout
and colour scheme are reproduced as closely as matplotlib allows.
=========================================================================
"""

from itertools import combinations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.interpolate import PchipInterpolator
    _HAS_PCHIP = True
except ImportError:                      # graceful fallback: linear ramp
    _HAS_PCHIP = False


# -----------------------------------------------------------------------
# Global font settings (Times New Roman with metric-compatible fallbacks)
# -----------------------------------------------------------------------
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.serif"] = ["Times New Roman", "Liberation Serif",
                              "TeX Gyre Termes", "DejaVu Serif"]
mpl.rcParams["mathtext.fontset"] = "stix"     # Times-like maths glyphs
mpl.rcParams["axes.unicode_minus"] = False


# -----------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------
BASELINE_MM = 130.16          # healthy geometry, healthy material

# Applied eta values at the M level (B level is zero for the four eta factors)
ETA_M = np.array([0.46, 0.25, 0.25, 1.655])   # [eta_D, eta_a, eta_T, eta_mu]
THETA_H = np.array([1.0, 2.0])                # [no hypertrophy, hypertrophy]

# Coded design in Yates order: -1 = B, +1 = M.
# Columns: run, x_D, x_a, x_T, x_mu, x_h, MM
raw = np.array([
    [1,  -1, -1, -1, -1, -1, 130.16],
    [2,  +1, -1, -1, -1, -1, 128.97],
    [3,  -1, +1, -1, -1, -1, 125.74],
    [4,  +1, +1, -1, -1, -1, 123.32],
    [5,  -1, -1, +1, -1, -1, 111.37],
    [6,  +1, -1, +1, -1, -1, 110.85],
    [7,  -1, +1, +1, -1, -1, 108.64],
    [8,  +1, +1, +1, -1, -1, 107.38],
    [9,  -1, -1, -1, +1, -1,  69.11],
    [10, +1, -1, -1, +1, -1,  67.91],
    [11, -1, +1, -1, +1, -1,  66.63],
    [12, +1, +1, -1, +1, -1,  66.42],
    [13, -1, -1, +1, +1, -1,  54.53],
    [14, +1, -1, +1, +1, -1,  54.19],
    [15, -1, +1, +1, +1, -1,  53.34],
    [16, +1, +1, +1, +1, -1,  53.15],
    [17, -1, -1, -1, -1, +1, 113.95],
    [18, +1, -1, -1, -1, +1, 114.44],
    [19, -1, +1, -1, -1, +1, 112.49],
    [20, +1, +1, -1, -1, +1, 111.04],
    [21, -1, -1, +1, -1, +1,  99.58],
    [22, +1, -1, +1, -1, +1,  99.31],
    [23, -1, +1, +1, -1, +1,  96.50],
    [24, +1, +1, +1, -1, +1,  95.97],
    [25, -1, -1, -1, +1, +1,  62.49],
    [26, +1, -1, -1, +1, +1,  62.32],
    [27, -1, +1, -1, +1, +1,  60.97],
    [28, +1, +1, -1, +1, +1,  60.84],
    [29, -1, -1, +1, +1, +1,  50.14],
    [30, +1, -1, +1, +1, +1,  49.94],
    [31, -1, +1, +1, +1, +1,  48.91],
    [32, +1, +1, +1, +1, +1,  48.80],
])

X = raw[:, 1:6]                # coded levels, +/-1
MM = raw[:, 6]                 # motility metric
n = raw.shape[0]

# Percentage change from the healthy baseline
PCT = (MM - BASELINE_MM) / BASELINE_MM * 100

# Log response
LOGMM = np.log(MM)

# Applied eta value for each run (0 at B, ETA_M at M)
ETA = (X[:, 0:4] + 1) / 2 * np.tile(ETA_M, (n, 1))

factors = ["eta_D", "eta_a", "eta_T", "eta_mu", "theta_h"]

# Maths bodies without the surrounding $ ... $, so that interaction labels
# can be built by concatenation later on.
math_bodies = [r"\eta_D", r"\eta_a", r"\eta_T", r"\eta_\mu",
               r"\vartheta^{\mathrm{h}}"]
labels = [f"${b}$" for b in math_bodies]

desc = ["Electrical diffusivity",
        "Excitation threshold",
        "Peak active stress",
        "Wall shear modulus",
        "Hypertrophy"]


def lvl(v):
    """Coded level -> 'B' (low) or 'M' (high)."""
    return "M" if v > 0 else "B"


# -----------------------------------------------------------------------
# Console: run table
# -----------------------------------------------------------------------
print("\n" + "=" * 78)
print(f"  Healthy baseline motility metric = {BASELINE_MM:.4f}")
print(f"  Levels: B (eta = 0, vartheta^h={THETA_H[0]:.1f})  |  "
      f"M (eta_D={ETA_M[0]:.3f}, eta_a={ETA_M[1]:.3f}, "
      f"eta_T={ETA_M[2]:.3f}, eta_mu={ETA_M[3]:.3f}, "
      f"vartheta^h={THETA_H[1]:.1f})")
print("=" * 78)
print(f"  {'Run':>3s}  {'eta_D':>6s} {'eta_a':>6s} {'eta_T':>6s} "
      f"{'eta_mu':>7s} {'theta_h':>8s}  {'MM':>9s}  {'pct change':>12s}")
print("-" * 78)

for i in range(n):
    print(f"  {int(raw[i, 0]):3d}  {lvl(X[i, 0]):>6s} {lvl(X[i, 1]):>6s} "
          f"{lvl(X[i, 2]):>6s} {lvl(X[i, 3]):>7s} {lvl(X[i, 4]):>8s}  "
          f"{MM[i]:9.4f}  {PCT[i]:+11.2f}%")

print("=" * 78)
imin = int(np.argmin(PCT))
imax = int(np.argmax(PCT))
print(f"  Min:  {PCT[imin]:+.2f}%  (Run {imin + 1:d})")
print(f"  Max:  {PCT[imax]:+.2f}%  (Run {imax + 1:d})")
print(f"  Mean: {PCT.mean():+.2f}%")
print("=" * 78 + "\n")


# -----------------------------------------------------------------------
# Regression: main effects + two-factor interactions
# -----------------------------------------------------------------------
pair_idx = np.array(list(combinations(range(5), 2)))   # 0-based factor pairs
n_pairs = pair_idx.shape[0]

inter = np.column_stack([X[:, i] * X[:, j] for i, j in pair_idx])

Xreg = np.column_stack([np.ones(n), X, inter])
beta = np.linalg.lstsq(Xreg, PCT, rcond=None)[0]

# Same model on the log scale, for the multiplicativity diagnostic
beta_log = np.linalg.lstsq(Xreg, LOGMM, rcond=None)[0]

coef_names = ["Intercept"] + factors
coef_names += [f"{factors[i]} x {factors[j]}" for i, j in pair_idx]

print("=" * 78)
print("REGRESSION MODEL  (coded inputs x_i in {-1,+1}; full effect = 2*beta)")
print("=" * 78)
print(f"{'Term':<24s} {'Coefficient':>14s} {'Full effect':>14s}")
print("-" * 78)

for k, b in enumerate(beta):
    if k == 0:
        print(f"{coef_names[k]:<24s} {b:+14.3f} {'-':>14s}")
    else:
        print(f"{coef_names[k]:<24s} {b:+14.3f} {2 * b:+14.3f}")

print("=" * 78)


# -----------------------------------------------------------------------
# Ranked effects
# -----------------------------------------------------------------------
term_names = coef_names[1:]
coeffs = beta[1:]
coeffs_log = beta_log[1:]
idx = np.argsort(-np.abs(coeffs), kind="stable")

print("\n" + "=" * 78)
print("RANKED EFFECTS")
print("=" * 78)
print(f"{'Rank':>4s} {'Term':<24s} {'Coeff (%)':>13s} "
      f"{'Effect (%)':>13s} {'log effect':>13s}")
print("-" * 78)

for k, j in enumerate(idx, start=1):
    print(f"{k:4d} {term_names[j]:<24s} {coeffs[j]:+13.3f} "
          f"{2 * coeffs[j]:+13.3f} {2 * coeffs_log[j]:+13.4f}")

print("=" * 78)


# -----------------------------------------------------------------------
# Main effects
# -----------------------------------------------------------------------
print("\n" + "=" * 78)
print("MAIN EFFECTS")
print("=" * 78)
print(f"{'Factor':<10s} {'Description':<24s} {'B mean':>10s} "
      f"{'M mean':>10s} {'Effect':>10s} {'MM ratio':>10s}")
print("-" * 78)

for k in range(5):
    lo_m = PCT[X[:, k] == -1].mean()
    hi_m = PCT[X[:, k] == +1].mean()
    ratio = MM[X[:, k] == +1].mean() / MM[X[:, k] == -1].mean()

    print(f"{factors[k]:<10s} {desc[k]:<24s} {lo_m:10.2f} {hi_m:10.2f} "
          f"{hi_m - lo_m:+10.2f} {ratio:10.4f}")

print("=" * 78)


# -----------------------------------------------------------------------
# Two-factor interactions
# -----------------------------------------------------------------------
print("\n" + "=" * 78)
print("TWO-FACTOR INTERACTIONS")
print("=" * 78)
print(f"{'Interaction':<26s} {'Coeff (%)':>14s} {'log coeff':>14s}")
print("-" * 78)

for p, (i, j) in enumerate(pair_idx):
    name = f"{factors[i]} x {factors[j]}"
    print(f"{name:<26s} {beta[6 + p]:+14.3f} {beta_log[6 + p]:+14.4f}")

print("=" * 78)


# -----------------------------------------------------------------------
# Multiplicativity check
# -----------------------------------------------------------------------
main_pct = np.mean(np.abs(beta[1:6]))
int_pct = np.mean(np.abs(beta[6:]))
main_log = np.mean(np.abs(beta_log[1:6]))
int_log = np.mean(np.abs(beta_log[6:]))

pred32 = MM[0]
for k in range(5):
    pred32 *= MM[X[:, k] == +1].mean() / MM[X[:, k] == -1].mean()

print("\n" + "=" * 78)
print("DIAGNOSTICS")
print("=" * 78)
print(f"  Interaction / main-effect ratio, pct scale : {int_pct / main_pct:.4f}")
print(f"  Interaction / main-effect ratio, log scale : {int_log / main_log:.4f}")

if int_log / main_log < 0.6 * int_pct / main_pct:
    print("  -> Interactions largely vanish on the log scale. They are "
          "scale-dependent")
    print("     rather than evidence of strong mechanistic synergy. Report "
          "the five")
    print("     mechanisms as acting independently and combining "
          "multiplicatively.")

print(f"\n  Multiplicative prediction of Run 32 from Run 1: {pred32:.4f}")
print(f"  Actual Run 32                                : {MM[31]:.4f}")
print("=" * 78 + "\n")


# =========================================================================
# Helper functions
# =========================================================================
def hex2rgb_local(hex_list):
    """Convert a list of '#rrggbb' strings to an (N, 3) array in [0, 1]."""
    return np.array([[int(h.lstrip("#")[i:i + 2], 16) / 255
                      for i in (0, 2, 4)] for h in hex_list])


def make_ramp(anchors, n_lev):
    """Smooth (pchip) interpolation of colour anchors onto n_lev levels."""
    xa = np.linspace(0, 1, anchors.shape[0])
    xq = np.linspace(0, 1, n_lev)

    if _HAS_PCHIP:
        cols = [PchipInterpolator(xa, anchors[:, c])(xq) for c in range(3)]
    else:
        cols = [np.interp(xq, xa, anchors[:, c]) for c in range(3)]

    return np.clip(np.column_stack(cols), 0, 1)


def ramp_color(frac, cmap):
    """Sample a colour from an (N, 3) ramp at fractional position frac."""
    frac = min(max(frac, 0.0), 1.0)
    return cmap[int(round(frac * (cmap.shape[0] - 1)))]


def set_figure_template(ax, font_size):
    """Common axis styling: serif ticks, boxed axes, square aspect, no grid."""
    ax.tick_params(labelsize=font_size, width=1.44,
                   direction="in", top=True, right=True)

    for spine in ax.spines.values():
        spine.set_linewidth(1.44)

    ax.grid(False)
    ax.set_box_aspect(1)


# =========================================================================
# COLOUR SCHEMES
# =========================================================================
RED_HEX = ["#f4e3d7", "#ffd5d5", "#ff8080", "#ff2a2a", "#aa0000"]
BLUE_HEX = ["#d7e3f4", "#d5e2ff", "#80a8ff", "#2a6aff", "#00368f"]

RED_ANCH = hex2rgb_local(RED_HEX)
BLUE_ANCH = hex2rgb_local(BLUE_HEX)

# Diverging map: dark blue -> light -> dark red
cmap_div = np.vstack([np.flipud(make_ramp(BLUE_ANCH, 128)),
                      make_ramp(RED_ANCH, 128)])

COL_B = RED_ANCH[1]
COL_M = RED_ANCH[2:4].mean(axis=0)
COL_PTS = np.array([0.67, 0.67, 0.67])
COL_LINE = np.array([0.07, 0.07, 0.07])

# Shared axis limits for the percentage plots
YLIM_PCT = (-80, 20)


# =========================================================================
# GLOBAL FIGURE FONT SIZES
# =========================================================================
# Increased throughout for improved readability

FS_TICK = 22          # axis tick numbers and B/M labels
FS_AXIS = 23          # x/y axis labels
FS_PANELTITLE = 22    # individual subplot titles
FS_TITLE = 28         # main figure titles
FS_LEGEND = 18        # legends


# -----------------------------------------------------------------------
# FIGURE 1: Diverging bar chart of regression coefficients
# LARGE-FONT PUBLICATION VERSION
# -----------------------------------------------------------------------
all_coeffs = beta[1:]

inter_labels_tex = [f"${math_bodies[i]} \\times {math_bodies[j]}$"
                    for i, j in pair_idx]
all_labels_tex = labels + inter_labels_tex

# Full ranking by magnitude
sidx_all = np.argsort(-np.abs(all_coeffs), kind="stable")
n_all = all_coeffs.size

# Only show coefficients above this magnitude
COEF_MIN = 0.5

keep = np.abs(all_coeffs[sidx_all]) >= COEF_MIN
seff = all_coeffs[sidx_all[keep]]
slabs = [all_labels_tex[k] for k in sidx_all[keep]]
n_show = seff.size
n_drop = n_all - n_show

print(f"Figure 1: showing {n_show} of {n_all} terms "
      f"(|beta| >= {COEF_MIN:.2f}%); {n_drop} omitted.")

# -----------------------------------------------------------------------
# Discrete colour scale
# -----------------------------------------------------------------------
CLIM_STEP = 5
CB_LABEL_STEP = 10

vabs = np.max(np.abs(all_coeffs))
CLIM = int(np.ceil(vabs / CLIM_STEP) * CLIM_STEP)
nBands = 2 * CLIM // CLIM_STEP

cmap_disc = np.array([ramp_color((k + 0.5) / nBands, cmap_div)
                      for k in range(nBands)])

x_lim = CLIM


def band_of(v):
    """Index of the discrete colour band containing coefficient v."""
    return int(min(max(np.floor((v + CLIM) / CLIM_STEP), 0), nBands - 1))


def col_fun(v):
    return cmap_disc[band_of(v)]


# -----------------------------------------------------------------------
# LARGE FONT SETTINGS - FIGURE 1
# -----------------------------------------------------------------------
FS_TICK_F1 = 36       # x-axis numbers + y-axis factor labels
FS_XLABEL_F1 = 37     # Regression coefficient (%)
FS_TITLE_F1 = 42      # main title
FS_ANNOT_F1 = 25      # negative / positive contribution text
                      # (sized so both annotations fit inside the axes)

FS_KEYTICK_F1 = 33    # colourbar numbers
FS_KEYTITLE_F1 = 34   # colourbar title

# -----------------------------------------------------------------------
# Large figure canvas
# -----------------------------------------------------------------------
FIG1_SIZE = (17.0, 9.0)     # 1700 x 900 px at 100 dpi
BAR_W = 0.42

fig1 = plt.figure(figsize=FIG1_SIZE, dpi=100, facecolor="w")

# Plot area is placed manually to leave room for the large y labels on the
# left and the discrete colour key on the right.
AX_POS = [0.155, 0.135, 0.585, 0.665]     # [left, bottom, width, height]
ax1 = fig1.add_axes(AX_POS)

# Zero line
ax1.axvline(0, color=(0.5, 0.5, 0.5), linewidth=1.8, zorder=1)

# Bars (rank 1 sits at the bottom, as in the MATLAB version)
for i, v in enumerate(seff, start=1):
    ax1.barh(i, v, height=BAR_W, color=col_fun(v),
             edgecolor="k", linewidth=1.8, zorder=2)

# Direction annotations
ax1.text(-x_lim * 0.97, n_show + 0.70,
         r"$\leftarrow$ Negative contribution to MM",
         fontsize=FS_ANNOT_F1, color=(0.25, 0.25, 0.25),
         ha="left", va="center")

ax1.text(x_lim * 0.97, n_show + 0.70,
         r"Positive contribution to MM $\rightarrow$",
         fontsize=FS_ANNOT_F1, color=(0.25, 0.25, 0.25),
         ha="right", va="center")

# Axes
ax1.set_yticks(np.arange(1, n_show + 1))
ax1.set_yticklabels(slabs)
ax1.set_xlim(-x_lim, x_lim)
ax1.set_ylim(0.30, n_show + 1.05)

ax1.tick_params(labelsize=FS_TICK_F1, width=1.8, length=8,
                direction="in", top=True, right=True)

for spine in ax1.spines.values():
    spine.set_linewidth(1.8)

ax1.set_xlabel("Regression coefficient (%)", fontsize=FS_XLABEL_F1)
ax1.grid(False)

# Title sits at the top of the canvas so that it clears the colour key
fig1.suptitle("Relative Contributions of Disease Parameters to Motility Loss",
              fontsize=FS_TITLE_F1, fontweight="bold", y=0.965)

# -----------------------------------------------------------------------
# Manual discrete colour key
# -----------------------------------------------------------------------
cbX = AX_POS[0] + AX_POS[2] + 0.035
cbY = AX_POS[1]
cbH = AX_POS[3]

bandH = cbH / nBands

# Make the bands approximately square in figure coordinates
cbW = bandH * FIG1_SIZE[1] / FIG1_SIZE[0]

cax = fig1.add_axes([cbX, cbY, cbW, cbH])
cax.set_xlim(0, 1)
cax.set_ylim(-CLIM, CLIM)

# Colour bands
for k in range(nBands):
    cax.axhspan(-CLIM + k * CLIM_STEP, -CLIM + (k + 1) * CLIM_STEP,
                facecolor=cmap_disc[k], edgecolor="none")

# Black separators between colour bands
for k in range(1, nBands):
    cax.axhline(-CLIM + k * CLIM_STEP, color="k", linewidth=1.0)

# Outer border of the colour key
for spine in cax.spines.values():
    spine.set_linewidth(1.2)
    spine.set_edgecolor("k")

# Tick marks every CLIM_STEP, numerical labels every CB_LABEL_STEP
ticks = np.arange(-CLIM, CLIM + CLIM_STEP, CLIM_STEP)
cax.set_xticks([])
cax.set_yticks(ticks)
cax.set_yticklabels([f"{t:+d}%" if t % CB_LABEL_STEP == 0 else ""
                     for t in ticks])
cax.yaxis.tick_right()
cax.tick_params(axis="y", labelsize=FS_KEYTICK_F1, width=1.1,
                length=6, direction="out")

# Colour-key title
cax.set_title("Regression\ncoefficient (%)", fontsize=FS_KEYTITLE_F1,
              fontweight="bold", pad=18)


# -----------------------------------------------------------------------
# FIGURE 2: Main effects (2 x 3)
# -----------------------------------------------------------------------
fig2 = plt.figure(figsize=(14.5, 9.0), dpi=100, facecolor="w")

for k in range(5):

    ax = fig2.add_subplot(2, 3, k + 1)

    lo_m = PCT[X[:, k] == -1].mean()
    hi_m = PCT[X[:, k] == +1].mean()

    ax.scatter(np.zeros(16), PCT[X[:, k] == -1],
               s=120, color=COL_PTS, alpha=0.85, zorder=2)

    ax.scatter(np.ones(16), PCT[X[:, k] == +1],
               s=120, color=COL_PTS, alpha=0.85, zorder=2)

    ax.plot([0, 1], [lo_m, hi_m], "o-", color=COL_LINE, linewidth=4.5,
            markersize=9, markerfacecolor=COL_LINE,
            markeredgecolor=COL_LINE, zorder=3)

    ax.axhline(0, linestyle=":", color=(0.55, 0.55, 0.55),
               linewidth=1.2, zorder=1)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["B", "M"])
    ax.set_xlim(-0.4, 1.4)
    ax.set_ylim(*YLIM_PCT)

    ax.set_title(f"{labels[k]}\nEffect = {hi_m - lo_m:+.1f}%",
                 fontsize=FS_PANELTITLE, fontweight="bold")

    ax.set_xlabel(labels[k], fontsize=FS_AXIS)

    if k % 3 == 0:
        ax.set_ylabel("% points", fontsize=FS_AXIS, rotation=0,
                      ha="right", va="center")

    set_figure_template(ax, FS_TICK)

fig2.suptitle("Main Effects on Motility Metric",
              fontsize=FS_TITLE, fontweight="bold")
fig2.tight_layout(rect=[0, 0, 1, 0.95], h_pad=3.0)


# -----------------------------------------------------------------------
# FIGURE 3: Interaction plots (2 x 5)
# -----------------------------------------------------------------------
LW_LINE = 4.0
LW_EDGE = 0.25
MS_PT = 200

# Enlarged canvas because the fonts are now larger
fig3 = plt.figure(figsize=(21.0, 9.5), dpi=100, facecolor="w")

for p, (ii, jj) in enumerate(pair_idx):

    ax = fig3.add_subplot(2, 5, p + 1)

    yB = [PCT[(X[:, jj] == -1) & (X[:, ii] == -1)].mean(),
          PCT[(X[:, jj] == -1) & (X[:, ii] == +1)].mean()]

    yM = [PCT[(X[:, jj] == +1) & (X[:, ii] == -1)].mean(),
          PCT[(X[:, jj] == +1) & (X[:, ii] == +1)].mean()]

    # Black underlay
    ax.plot([0, 1], yB, "-", color="k", linewidth=LW_LINE + 2 * LW_EDGE,
            zorder=2)
    ax.plot([0, 1], yM, "-", color="k", linewidth=LW_LINE + 2 * LW_EDGE,
            zorder=2)

    # Coloured lines
    ax.plot([0, 1], yB, "-", color=COL_B, linewidth=LW_LINE,
            label=f"{labels[jj]} = B", zorder=3)
    ax.plot([0, 1], yM, "-", color=COL_M, linewidth=LW_LINE,
            label=f"{labels[jj]} = M", zorder=3)

    # Markers
    ax.scatter([0, 1], yB, s=MS_PT, color=COL_B, edgecolor="k",
               linewidth=LW_EDGE, zorder=4)
    ax.scatter([0, 1], yM, s=MS_PT, color=COL_M, edgecolor="k",
               linewidth=LW_EDGE, zorder=4)

    ax.axhline(0, linestyle=":", color=(0.55, 0.55, 0.55),
               linewidth=1.2, zorder=1)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["B", "M"])
    ax.set_xlim(-0.3, 1.3)
    ax.set_ylim(*YLIM_PCT)

    ax.set_title(f"${math_bodies[ii]} \\times {math_bodies[jj]}$\n"
                 f"coeff = {beta[6 + p]:+.2f}%",
                 fontsize=FS_PANELTITLE, fontweight="bold")

    ax.set_xlabel(labels[ii], fontsize=FS_AXIS)

    if p % 5 == 0:
        ax.set_ylabel("% points", fontsize=FS_AXIS, rotation=0,
                      ha="right", va="center")

    ax.legend(fontsize=FS_LEGEND, loc="lower left", frameon=False)

    set_figure_template(ax, FS_TICK)

fig3.suptitle("Two-Factor Interaction Plots",
              fontsize=FS_TITLE, fontweight="bold")
fig3.tight_layout(rect=[0, 0, 1, 0.94], h_pad=3.5)

plt.show()
