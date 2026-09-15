"""
tune_fiber_stiffness.py

Standalone SEQUENTIAL calibration of an anisotropic hyperelastic model
against two uniaxial stress-stretch datasets.

Place these three files in the same folder:
    tune_fiber_stiffness.py
    circum_stretch_stress.txt
    long_stretch_stress.txt


WORKFLOW
--------
Stage 1   Fit KF_CIRC against the circumferential dataset alone.
          KF_LONG is held at KF_LONG_PLACEHOLDER for this stage.

Stage 2   Fit KF_LONG against the longitudinal dataset alone,
          with KF_CIRC held fixed at the Stage 1 result.

Stage 3   Consistency check. Stage 1 is repeated using the Stage 2
          value of KF_LONG. If KF_CIRC is unchanged, the two stages do
          not interact and the sequential fit is exact rather than an
          approximation to a joint fit.


MATERIAL MODEL
--------------
Strain energy (compressible neo-Hookean matrix plus two tension-only
fibre families with a cubic engagement law):

    W = (mu0/2)(I1 - 3) - mu0*ln(J) + (kappa/2)(ln J)^2
      + (k_c/2) <I4c - 1>^3
      + (k_l/2) <I4l - 1>^3

with kappa = KAPPA_FACTOR * mu0 and <x> = max(x, 0).

Cauchy stress:

    sigma = (mu0/J)(B - I) + (kappa*ln(J)/J) I
          + (3*k_c/J) I4c <I4c - 1>^2 (e_c (x) e_c)
          + (3*k_l/J) I4l <I4l - 1>^2 (e_l (x) e_l)

Fibre 2 (e2) is circumferential, fibre 3 (e3) is longitudinal.
The model is exactly stress-free at F = I.


TESTS REPRODUCED
----------------
Circumferential   imposed lambda_2, solve sigma_11 = sigma_33 = 0,
                  fit sigma_22
Longitudinal      imposed lambda_3, solve sigma_11 = sigma_22 = 0,
                  fit sigma_33


OUTPUTS
-------
    best_fit_parameters.txt
    sequential_fit_metrics.txt
    circum_python_fit.txt   circum_python_fit_metrics.txt   circum_python_fit.png
    long_python_fit.txt     long_python_fit_metrics.txt     long_python_fit.png

Requires: numpy  scipy  matplotlib
"""

from pathlib import Path
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares


# =====================================================================
# USER SETTINGS
# =====================================================================

CIRCUM_FILE = "circum_stretch_stress.txt"
LONG_FILE = "long_stretch_stress.txt"

# Starting guesses for the two tuned parameters.
KF_CIRC_INITIAL = 210.0
KF_LONG_INITIAL = 50.0

# Value used for KF_LONG during Stage 1, before it has been fitted.
# Stage 3 verifies whether this choice mattered.
KF_LONG_PLACEHOLDER = KF_LONG_INITIAL

# Matrix parameters are held fixed throughout.
MU0_FIXED = 0.64
KAPPA_FACTOR = 100.0

# Optional baseline subtraction, applied to each dataset separately.
#   False = fit the stresses exactly as supplied
#   True  = subtract each dataset's first stress value
#
# The model is stress-free at lambda = 1, so any non-zero stress at the
# first data point is an offset that no choice of k_c or k_l can match.
BASELINE_CORRECT = False

# Broad positive bounds for the tuned stiffnesses.
KF_MIN = 1.0e-8
KF_MAX = 1.0e8

# Admissible range for the solved lateral stretches.
LATERAL_MIN = 0.02
LATERAL_MAX = 5.0

# A lateral-equilibrium solve is accepted when
#   max|sigma_lateral| / max(|sigma_axial|, matrix stiffness) < this.
LATERAL_TOLERANCE = 1.0e-8

# Stage 3 flags a coupling problem if KF_CIRC moves by more than this
# relative amount when Stage 1 is repeated with the fitted KF_LONG.
CONSISTENCY_TOLERANCE = 1.0e-6


# =====================================================================
# DERIVED CONSTANTS
# =====================================================================

KAPPA_B = KAPPA_FACTOR * MU0_FIXED

# Stiffness governing the lateral equilibrium equations, used to
# normalise the lateral residuals.
#
# This must NOT be scaled by k_c or k_l: those can reach 1e8 under the
# bounds above, and dividing by them drives the normalised residual
# gradient below gtol, so the solver stops on its initial guess and
# silently returns a state that is not in equilibrium.
MATRIX_SCALE = MU0_FIXED + KAPPA_B

# Per-direction configuration, replacing branched 'circum'/'long' logic.
#
#   axial    index of the imposed-stretch axis
#   lateral  indices of the two traction-free axes
#   labels   column headings for the two lateral stretches
DIRECTIONS = {
    "circum": {
        "axial": 1,
        "lateral": (0, 2),
        "labels": ("lambda_x", "lambda_z"),
        "stress_symbol": r"$\sigma_{22}$",
        "title": "Circumferential Passive Response",
        "data_file": "circum_python_fit.txt",
        "metric_file": "circum_python_fit_metrics.txt",
        "plot_file": "circum_python_fit.png",
    },
    "long": {
        "axial": 2,
        "lateral": (0, 1),
        "labels": ("lambda_x", "lambda_y"),
        "stress_symbol": r"$\sigma_{33}$",
        "title": "Longitudinal Passive Response",
        "data_file": "long_python_fit.txt",
        "metric_file": "long_python_fit_metrics.txt",
        "plot_file": "long_python_fit.png",
    },
}


# =====================================================================
# FILE LOCATION
# =====================================================================

try:
    SCRIPT_DIRECTORY = Path(__file__).resolve().parent
except NameError:
    # Interactive session: fall back to the working directory.
    SCRIPT_DIRECTORY = Path.cwd()

CIRCUM_PATH = SCRIPT_DIRECTORY / CIRCUM_FILE
LONG_PATH = SCRIPT_DIRECTORY / LONG_FILE


# =====================================================================
# READ EXPERIMENTAL DATA
# =====================================================================

def read_experimental_data(path):
    """
    Read a two-column text file of  stretch  stress  and return both
    columns sorted by ascending stretch.
    """

    if not path.exists():
        raise FileNotFoundError(
            f"\nCould not find experimental file:\n\n  {path}\n\n"
            "Put it in the same folder as this script."
        )

    data = np.loadtxt(path, dtype=float)

    # A single-row file loads as 1-D.
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 2:
        raise ValueError(
            f"{path.name} must contain at least two columns: stretch stress"
        )

    stretch = data[:, 0].astype(float)
    stress = data[:, 1].astype(float)

    if stretch.size < 2:
        raise ValueError(f"{path.name} needs at least two data points.")

    if np.any(stretch <= 0.0):
        raise ValueError(f"{path.name} contains a non-positive stretch.")

    order = np.argsort(stretch)
    stretch = stretch[order]
    stress = stress[order]

    if BASELINE_CORRECT:
        baseline = stress[0]
        stress = stress - baseline
        print(f"Baseline-corrected {path.name} by subtracting {baseline:.8f}")

    return stretch, stress


# =====================================================================
# CAUCHY STRESS
# =====================================================================

def cauchy_stress(F, kf_circ, kf_long):
    """
    Cauchy stress for a given deformation gradient.

    Parameters
    ----------
    F        3x3 deformation gradient
    kf_circ  circumferential fibre stiffness k_c
    kf_long  longitudinal fibre stiffness k_l

    Returns
    -------
    3x3 Cauchy stress tensor.
    """

    F = np.asarray(F, dtype=float).reshape(3, 3)

    J = float(np.linalg.det(F))
    if J <= 1.0e-12:
        J = 1.0e-12

    logJ = math.log(J)

    # Left Cauchy-Green tensor.
    B = F @ F.T

    # Matrix contribution. At F = I this gives mu0*I - mu0*I = 0, so the
    # reference state is stress-free.
    p_scalar = KAPPA_B * logJ - MU0_FIXED
    sigma = (MU0_FIXED / J) * B + (p_scalar / J) * np.eye(3)

    # Fibre contributions. e2 = circumferential, e3 = longitudinal.
    fibres = (
        (np.array([0.0, 1.0, 0.0]), kf_circ),
        (np.array([0.0, 0.0, 1.0]), kf_long),
    )

    for direction_ref, kf in fibres:

        # Push the fibre forward and extract its stretch.
        a = F @ direction_ref
        lam = max(float(np.linalg.norm(a)), 1.0e-12)
        e = a / lam

        i4 = lam ** 2

        # Tension-only: a compressed fibre carries no stress.
        bracket = max(i4 - 1.0, 0.0)

        if bracket > 0.0:
            sigma = sigma + (
                (3.0 * kf / J) * i4 * bracket ** 2 * np.outer(e, e)
            )

    return sigma


def fibre_stretches(F):
    """Return (lambda_circ, lambda_long) for a deformation gradient."""

    F = np.asarray(F, dtype=float).reshape(3, 3)

    return (
        float(np.linalg.norm(F @ np.array([0.0, 1.0, 0.0]))),
        float(np.linalg.norm(F @ np.array([0.0, 0.0, 1.0]))),
    )


# =====================================================================
# SOLVE ONE UNIAXIAL STATE
# =====================================================================

def _build_F(applied_stretch, lateral, config):
    """Assemble the diagonal deformation gradient for one test."""

    stretches = np.empty(3, dtype=float)
    stretches[config["axial"]] = applied_stretch
    stretches[config["lateral"][0]] = lateral[0]
    stretches[config["lateral"][1]] = lateral[1]

    return np.diag(stretches)


def solve_uniaxial_state(
        applied_stretch,
        direction,
        kf_circ,
        kf_long,
        lateral_guess=None):
    """
    Solve for the lateral stretches that make both traction-free faces
    stress-free at a prescribed axial stretch.

    The unknowns are log(lateral stretches), which keeps them positive.

    Returns
    -------
    axial_stress     stress along the imposed direction
    lateral_1        first solved lateral stretch
    lateral_2        second solved lateral stretch
    J                det(F)
    equilibrium_err  max|sigma_lateral| normalised by the local stress
                     scale; should be ~1e-12 or smaller
    """

    if direction not in DIRECTIONS:
        raise ValueError(
            f"direction must be one of {sorted(DIRECTIONS)}, got {direction!r}"
        )

    config = DIRECTIONS[direction]
    i_lat_a, i_lat_b = config["lateral"]
    i_axial = config["axial"]

    lam = float(applied_stretch)
    if lam <= 0.0:
        raise ValueError("Stretch must be positive.")

    # Incompressible estimate, used for the first point and for retries.
    cold_guess = np.full(2, lam ** (-0.5), dtype=float)

    if lateral_guess is None:
        lateral_guess = cold_guess
    else:
        lateral_guess = np.asarray(lateral_guess, dtype=float)

    def residual(log_lateral):
        """Normalised [sigma_lat_a, sigma_lat_b]."""
        lateral = np.exp(log_lateral)
        sigma = cauchy_stress(_build_F(lam, lateral, config), kf_circ, kf_long)
        return np.array([sigma[i_lat_a, i_lat_a],
                         sigma[i_lat_b, i_lat_b]]) / MATRIX_SCALE

    def run(guess):
        """Solve from one starting guess and report the true error."""

        # Keep the start strictly inside the bounds, otherwise
        # least_squares rejects it.
        x0 = np.log(np.clip(guess, LATERAL_MIN * 1.5, LATERAL_MAX / 1.5))

        solution = least_squares(
            residual,
            x0,
            bounds=(
                np.log([LATERAL_MIN, LATERAL_MIN]),
                np.log([LATERAL_MAX, LATERAL_MAX]),
            ),
            xtol=1.0e-14,
            ftol=1.0e-14,
            gtol=1.0e-14,
            max_nfev=300,
        )

        lateral = np.exp(solution.x)
        F = _build_F(lam, lateral, config)
        sigma = cauchy_stress(F, kf_circ, kf_long)

        lateral_stress = max(
            abs(sigma[i_lat_a, i_lat_a]),
            abs(sigma[i_lat_b, i_lat_b]),
        )

        # Measure the error against the stress actually present, so the
        # check means the same thing at every load level.
        err = lateral_stress / max(abs(sigma[i_axial, i_axial]), MATRIX_SCALE)

        return lateral, sigma, F, err

    lateral, sigma, F, err = run(lateral_guess)

    # A warm start inherited from a previous point can stall. Retry cold
    # before accepting a state that is not in equilibrium.
    if err > LATERAL_TOLERANCE and not np.allclose(lateral_guess, cold_guess):
        retry = run(cold_guess)
        if retry[3] < err:
            lateral, sigma, F, err = retry

    return (
        float(sigma[i_axial, i_axial]),
        float(lateral[0]),
        float(lateral[1]),
        float(np.linalg.det(F)),
        float(err),
    )


# =====================================================================
# CALCULATE COMPLETE CURVE
# =====================================================================

def calculate_model_curve(stretches, direction, kf_circ, kf_long):
    """
    Sweep a stretch list and return the full uniaxial response.

    Each point warm-starts from the previous converged lateral stretches,
    which makes the sweep cheaper and more stable.

    Returns
    -------
    stresses, lateral_1, lateral_2, J_values, worst_equilibrium_error
    """

    stretches = np.asarray(stretches, dtype=float)
    n = len(stretches)

    stresses = np.zeros(n, dtype=float)
    lateral_1 = np.zeros(n, dtype=float)
    lateral_2 = np.zeros(n, dtype=float)
    J_values = np.zeros(n, dtype=float)

    worst_error = 0.0
    previous_lateral = None

    for i, stretch in enumerate(stretches):

        axial_stress, lat1, lat2, J, err = solve_uniaxial_state(
            stretch, direction, kf_circ, kf_long, previous_lateral
        )

        stresses[i] = axial_stress
        lateral_1[i] = lat1
        lateral_2[i] = lat2
        J_values[i] = J

        worst_error = max(worst_error, err)
        previous_lateral = np.array([lat1, lat2], dtype=float)

    return stresses, lateral_1, lateral_2, J_values, worst_error


# =====================================================================
# METRICS
# =====================================================================

def calculate_metrics(experimental, predicted):
    """Return (RMSE, NRMSE, R^2). NRMSE is scaled by the peak |stress|."""

    experimental = np.asarray(experimental, dtype=float)
    predicted = np.asarray(predicted, dtype=float)

    rmse = math.sqrt(np.mean((predicted - experimental) ** 2))
    nrmse = rmse / max(float(np.max(np.abs(experimental))), 1.0)

    ss_res = np.sum((experimental - predicted) ** 2)
    ss_tot = np.sum((experimental - np.mean(experimental)) ** 2)

    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")

    return rmse, nrmse, r_squared


# =====================================================================
# SINGLE-PARAMETER FIT (ONE STAGE)
# =====================================================================

def fit_one_direction(
        stretch, stress, direction,
        initial_kf, other_kf, verbose=True):
    """
    Tune the fibre stiffness belonging to one direction.

    'circum' tunes KF_CIRC with KF_LONG held at other_kf.
    'long'   tunes KF_LONG with KF_CIRC held at other_kf.

    Optimisation is over log(k), which enforces positivity and makes the
    step size scale-free across the very wide bounds.

    Returns
    -------
    kf_best, result, evaluations
    """

    if direction not in DIRECTIONS:
        raise ValueError(
            f"direction must be one of {sorted(DIRECTIONS)}, got {direction!r}"
        )

    tuned_name = "KF_CIRC" if direction == "circum" else "KF_LONG"
    fixed_name = "KF_LONG" if direction == "circum" else "KF_CIRC"

    # Residuals are normalised by the peak measured stress so the
    # convergence tolerances mean the same thing for both datasets.
    scale = max(float(np.max(np.abs(stress))), 1.0)

    evaluations = [0]

    def objective(theta):

        kf = float(np.exp(theta[0]))
        kf_circ, kf_long = (
            (kf, other_kf) if direction == "circum" else (other_kf, kf)
        )

        model, _, _, _, _ = calculate_model_curve(
            stretch, direction, kf_circ, kf_long
        )

        evaluations[0] += 1

        if verbose and (evaluations[0] == 1 or evaluations[0] % 5 == 0):
            rmse = math.sqrt(np.mean((model - stress) ** 2))
            print(
                f"  Evaluation {evaluations[0]:4d}: "
                f"{tuned_name}={kf:14.6f}   RMSE={rmse:14.6f}"
            )

        return (model - stress) / scale

    result = least_squares(
        objective,
        np.log([initial_kf]),
        bounds=(np.log([KF_MIN]), np.log([KF_MAX])),
        xtol=1.0e-12,
        ftol=1.0e-12,
        gtol=1.0e-12,
        max_nfev=500,
        verbose=0,
    )

    return float(np.exp(result.x[0])), result, evaluations[0]


def off_axis_fibre_activates(stretch, lateral_1, lateral_2, direction):
    """
    Report whether the fibre family NOT aligned with the loading axis
    ever goes into tension.

    If it never does, that family contributes no stress to this test, so
    its stiffness has no effect on the fitted value.
    """

    config = DIRECTIONS[direction]

    for i in range(len(stretch)):
        F = _build_F(stretch[i], (lateral_1[i], lateral_2[i]), config)
        lam_c, lam_l = fibre_stretches(F)
        off_axis = lam_l if direction == "circum" else lam_c
        if off_axis > 1.0:
            return True

    return False


# =====================================================================
# WRITE ONE DIRECTION'S RESULTS
# =====================================================================

def save_direction_results(
        direction, stretch, stress_exp, stress_fit,
        lateral_1, lateral_2, J):
    """Write the per-point table and metrics file for one direction."""

    config = DIRECTIONS[direction]
    label_a, label_b = config["labels"]

    rmse, nrmse, r_squared = calculate_metrics(stress_exp, stress_fit)

    with open(SCRIPT_DIRECTORY / config["data_file"], "w") as f:
        f.write(f"# lambda sigma_exp sigma_fit {label_a} {label_b} J\n")
        for i in range(len(stretch)):
            f.write(
                f"{stretch[i]:.8f} {stress_exp[i]:.8f} {stress_fit[i]:.8f} "
                f"{lateral_1[i]:.8f} {lateral_2[i]:.8f} {J[i]:.8f}\n"
            )

    with open(SCRIPT_DIRECTORY / config["metric_file"], "w") as f:
        f.write(f"RMSE {rmse:.12f}\n")
        f.write(f"NRMSE {nrmse:.12f}\n")
        f.write(f"R2 {r_squared:.12f}\n")

    return rmse, nrmse, r_squared


# =====================================================================
# PLOT ONE DIRECTION
# =====================================================================

def make_plot(direction, stretch, stress_exp, stress_fit, show=True):
    """Save a data-versus-fit figure for one direction."""

    config = DIRECTIONS[direction]

    fig, ax = plt.subplots(figsize=(9, 7))

    ax.plot(
        stretch, stress_exp, "o",
        markersize=8,
        markerfacecolor="0.75",
        markeredgecolor="black",
        markeredgewidth=1.3,
        linestyle="None",
        label="Experimental Data",
    )

    ax.plot(
        stretch, stress_fit, "-",
        linewidth=3.5,
        color="black",
        label="Fitted Model",
    )

    ax.set_xlabel(r"Stretch $\lambda$", fontsize=24)
    ax.set_ylabel(f"Cauchy Stress {config['stress_symbol']}", fontsize=24)
    ax.set_title(config["title"], fontsize=24)

    ax.legend(fontsize=15, loc="upper left")
    ax.grid(False)
    ax.tick_params(labelsize=17, width=1.5, length=6)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    ax.set_xlim(float(np.min(stretch)), float(np.max(stretch)))

    fig.tight_layout()

    plot_path = SCRIPT_DIRECTORY / config["plot_file"]
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"Saved plot: {plot_path}")

    if show:
        plt.show()

    # Release the figure so repeated runs do not accumulate open figures.
    plt.close(fig)


# =====================================================================
# MAIN
# =====================================================================

def main(show_plots=True):

    print("")
    print("============================================================")
    print("SEQUENTIAL CALIBRATION:  CIRCUMFERENTIAL, THEN LONGITUDINAL")
    print("============================================================")
    print("")
    print(f"Script directory     : {SCRIPT_DIRECTORY}")
    print(f"Circumferential data : {CIRCUM_PATH}")
    print(f"Longitudinal data    : {LONG_PATH}")

    stretch_c, stress_c = read_experimental_data(CIRCUM_PATH)
    stretch_l, stress_l = read_experimental_data(LONG_PATH)

    print("")
    print(f"Circumferential points : {len(stretch_c)}")
    print(f"Longitudinal points    : {len(stretch_l)}")
    print("")
    print("Fixed matrix parameters:")
    print(f"  MU0     = {MU0_FIXED:.8f}")
    print(f"  KAPPA_B = {KAPPA_B:.8f}")

    # The model is stress-free at lambda = 1, so a non-zero first data
    # point leaves a residual that tuning cannot remove.
    if not BASELINE_CORRECT:
        for name, stretch, stress in (
            ("circumferential", stretch_c, stress_c),
            ("longitudinal", stretch_l, stress_l),
        ):
            if abs(stress[0]) > 1.0e-10:
                print(
                    f"\nNOTE: {name} stress at lambda = {stretch[0]:.4f} is "
                    f"{stress[0]:.6f}, but the model is stress-free at F = I."
                    "\n      This offset cannot be fitted. Consider "
                    "BASELINE_CORRECT = True."
                )

    # ================================================================
    # STAGE 1 - CIRCUMFERENTIAL
    # ================================================================

    print("")
    print("============================================================")
    print("STAGE 1 of 2  -  CIRCUMFERENTIAL")
    print("============================================================")
    print("")
    print(f"  Tuning  KF_CIRC   from initial guess {KF_CIRC_INITIAL:.6f}")
    print(f"  Holding KF_LONG   at placeholder     {KF_LONG_PLACEHOLDER:.6f}")
    print("")

    kf_circ_best, result_c, evals_c = fit_one_direction(
        stretch_c, stress_c, "circum",
        initial_kf=KF_CIRC_INITIAL,
        other_kf=KF_LONG_PLACEHOLDER,
    )

    fit_c, lat1_c, lat2_c, J_c, err_c = calculate_model_curve(
        stretch_c, "circum", kf_circ_best, KF_LONG_PLACEHOLDER
    )

    rmse_c, nrmse_c, r2_c = save_direction_results(
        "circum", stretch_c, stress_c, fit_c, lat1_c, lat2_c, J_c
    )

    long_fibre_active = off_axis_fibre_activates(
        stretch_c, lat1_c, lat2_c, "circum"
    )

    print("")
    print("  ------------------------------------------------------")
    print(f"  STAGE 1 RESULT:   KF_CIRC = {kf_circ_best:.10f}")
    print("  ------------------------------------------------------")
    print("")
    print(f"  RMSE  = {rmse_c:.10f}")
    print(f"  NRMSE = {nrmse_c:.10f}")
    print(f"  R^2   = {r2_c:.10f}")
    print("")
    print(f"  Optimiser success    = {result_c.success}")
    print(f"  Residual evaluations = {evals_c}")
    print(f"  Lateral equilibrium  = {err_c:.3e}")
    print(f"  J range              = {J_c.min():.6f} to {J_c.max():.6f}")
    print(
        "  Longitudinal fibre   = "
        + (
            "IN TENSION somewhere, so the KF_LONG placeholder "
            "influenced this result"
            if long_fibre_active
            else "in compression throughout, so KF_LONG had no "
                 "influence on this result"
        )
    )

    # ================================================================
    # STAGE 2 - LONGITUDINAL
    # ================================================================

    print("")
    print("============================================================")
    print("STAGE 2 of 2  -  LONGITUDINAL")
    print("============================================================")
    print("")
    print(f"  Tuning  KF_LONG   from initial guess {KF_LONG_INITIAL:.6f}")
    print(f"  Holding KF_CIRC   at Stage 1 result  {kf_circ_best:.6f}")
    print("")

    kf_long_best, result_l, evals_l = fit_one_direction(
        stretch_l, stress_l, "long",
        initial_kf=KF_LONG_INITIAL,
        other_kf=kf_circ_best,
    )

    fit_l, lat1_l, lat2_l, J_l, err_l = calculate_model_curve(
        stretch_l, "long", kf_circ_best, kf_long_best
    )

    rmse_l, nrmse_l, r2_l = save_direction_results(
        "long", stretch_l, stress_l, fit_l, lat1_l, lat2_l, J_l
    )

    circ_fibre_active = off_axis_fibre_activates(
        stretch_l, lat1_l, lat2_l, "long"
    )

    print("")
    print("  ------------------------------------------------------")
    print(f"  STAGE 2 RESULT:   KF_LONG = {kf_long_best:.10f}")
    print("  ------------------------------------------------------")
    print("")
    print(f"  RMSE  = {rmse_l:.10f}")
    print(f"  NRMSE = {nrmse_l:.10f}")
    print(f"  R^2   = {r2_l:.10f}")
    print("")
    print(f"  Optimiser success    = {result_l.success}")
    print(f"  Residual evaluations = {evals_l}")
    print(f"  Lateral equilibrium  = {err_l:.3e}")
    print(f"  J range              = {J_l.min():.6f} to {J_l.max():.6f}")
    print(
        "  Circumferential fibre = "
        + (
            "IN TENSION somewhere, so KF_CIRC influenced this result"
            if circ_fibre_active
            else "in compression throughout, so KF_CIRC had no "
                 "influence on this result"
        )
    )

    # ================================================================
    # STAGE 3 - CONSISTENCY CHECK
    # ================================================================
    #
    # Stage 1 had to assume a value for KF_LONG. Repeating it with the
    # fitted KF_LONG shows whether that assumption mattered. No change
    # means the two stages are independent and running them in sequence
    # gives exactly the same answer as fitting them together.

    print("")
    print("============================================================")
    print("STAGE 3  -  CONSISTENCY CHECK")
    print("============================================================")
    print("")
    print(f"  Refitting KF_CIRC with KF_LONG = {kf_long_best:.6f}")

    kf_circ_recheck, _, _ = fit_one_direction(
        stretch_c, stress_c, "circum",
        initial_kf=kf_circ_best,
        other_kf=kf_long_best,
        verbose=False,
    )

    drift = abs(kf_circ_recheck - kf_circ_best) / max(abs(kf_circ_best), 1e-30)

    print("")
    print(f"  KF_CIRC from Stage 1 : {kf_circ_best:.10f}")
    print(f"  KF_CIRC on refit     : {kf_circ_recheck:.10f}")
    print(f"  Relative change      : {drift:.3e}")
    print("")

    if drift <= CONSISTENCY_TOLERANCE:
        print("  The stages are independent. The sequential fit is exact,")
        print("  not an approximation, and a second pass is unnecessary.")
    else:
        print("  The stages interact: KF_CIRC moved once KF_LONG was known.")
        print("  Re-run Stage 1 and Stage 2 until the values settle, or fit")
        print("  both parameters together instead.")

    # ================================================================
    # SAVE PARAMETERS AND METRICS
    # ================================================================

    with open(SCRIPT_DIRECTORY / "best_fit_parameters.txt", "w") as f:
        f.write(f"KF_CIRC {kf_circ_best:.12f}\n")
        f.write(f"KF_LONG {kf_long_best:.12f}\n")
        f.write(f"MU0 {MU0_FIXED:.12f}\n")
        f.write(f"KAPPA_B {KAPPA_B:.12f}\n")
        f.write("FIT_MODE sequential_circum_then_long\n")
        f.write(f"STAGE1_SUCCESS {result_c.success}\n")
        f.write(f"STAGE2_SUCCESS {result_l.success}\n")
        f.write(f"STAGE1_EVALUATIONS {evals_c}\n")
        f.write(f"STAGE2_EVALUATIONS {evals_l}\n")
        f.write(f"CONSISTENCY_DRIFT {drift:.12e}\n")

    with open(SCRIPT_DIRECTORY / "sequential_fit_metrics.txt", "w") as f:
        f.write(f"CIRCUM_RMSE {rmse_c:.12f}\n")
        f.write(f"CIRCUM_NRMSE {nrmse_c:.12f}\n")
        f.write(f"CIRCUM_R2 {r2_c:.12f}\n")
        f.write(f"LONG_RMSE {rmse_l:.12f}\n")
        f.write(f"LONG_NRMSE {nrmse_l:.12f}\n")
        f.write(f"LONG_R2 {r2_l:.12f}\n")

    # ================================================================
    # SUMMARY
    # ================================================================

    print("")
    print("============================================================")
    print("FINAL PARAMETERS")
    print("============================================================")
    print("")
    print(f"  KF_CIRC = {kf_circ_best:.10f}   (Stage 1, circumferential)")
    print(f"  KF_LONG = {kf_long_best:.10f}   (Stage 2, longitudinal)")
    print(f"  MU0     = {MU0_FIXED:.10f}   (fixed)")
    print(f"  KAPPA_B = {KAPPA_B:.10f}   (fixed)")
    print("")
    print(f"  Circumferential   R^2 = {r2_c:.10f}   RMSE = {rmse_c:.6f}")
    print(f"  Longitudinal      R^2 = {r2_l:.10f}   RMSE = {rmse_l:.6f}")
    print("")
    print("Output files:")
    for name in (
        "best_fit_parameters.txt",
        "sequential_fit_metrics.txt",
        "circum_python_fit.txt",
        "circum_python_fit_metrics.txt",
        "circum_python_fit.png",
        "long_python_fit.txt",
        "long_python_fit_metrics.txt",
        "long_python_fit.png",
    ):
        print(f"  {name}")

    make_plot("circum", stretch_c, stress_c, fit_c, show=show_plots)
    make_plot("long", stretch_l, stress_l, fit_l, show=show_plots)


# =====================================================================
# RUN
# =====================================================================

if __name__ == "__main__":
    main()
