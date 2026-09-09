"""
extract_and_plot_stress_strain.py

Post-processor for single_elem_circum.inp / single_elem_long.inp.
Three modes, run in this order:

Step 1 (run with Abaqus python, has odbAccess but usually no matplotlib):
    abaqus python extract_and_plot_stress_strain.py --extract single_elem_circum circum

    Opens single_elem_circum.odb, pulls the history outputs for the loaded
    face (U2/RF2 for 'circum', U3/RF3 for 'long'), converts them to
    engineering stretch/stress, and writes a plain text file (columns:
    stretch, stress), e.g. circum_job_response.txt.

Step 2a - quick look (run with a normal python that has matplotlib):
    python extract_and_plot_stress_strain.py --plot circum_job_response.txt

    Plots the raw simulated curve only, saved as stress_strain.png.
    (If matplotlib isn't installed: pip install matplotlib)

Step 2b - publication-style figure with experimental overlay + fit metrics,
matching the style/content of the old score_circumi_wrapper.m /
score_longi_wrapper2.m MATLAB scripts (no MATLAB or Abaqus-Intel setup
needed):
    python extract_and_plot_stress_strain.py --score circum_job_response.txt circum_stretch_stress.txt circum
    python extract_and_plot_stress_strain.py --score long_job_response.txt  long_stretch_stress.txt  long

    This:
      - loads the experimental (stretch, stress) data
      - resamples it onto a 0.05-stretch-increment grid (linear, extrapolated)
      - interpolates the simulated response onto the same grid
      - computes RMSE / NRMSE on that grid
      - writes passive_circum_export.txt / passive_longi_export.txt
        (columns: lambda, sigma_exp, sigma_fit)
      - writes circumi_fit_metrics.txt / longi_fit_metrics.txt
      - saves a large, publication-style PNG (gray filled circles for
        Data, thick black line for Fitted Model, big serif/LaTeX-style
        labels, square box, no grid) as circum_publication_fit.png /
        long_publication_fit.png

Everything lives in one file so there is only one script to keep track of;
each mode only needs the libraries available in its own interpreter.
"""

import sys


def extract(jobname, mode):
    from odbAccess import openOdb

    odb = openOdb(jobname + '.odb')
    step = odb.steps[list(odb.steps.keys())[-1]]

    history_region = None
    for key in step.historyRegions.keys():
        hr = step.historyRegions[key]
        names = [k.upper() for k in hr.historyOutputs.keys()]
        if mode == 'circum':
            if 'U2' in names and 'RF2' in names:
                history_region = hr
                break
        elif mode == 'long':
            if 'U3' in names and 'RF3' in names:
                history_region = hr
                break

    if history_region is None:
        raise RuntimeError('Could not find appropriate history region.')

    u_key = None
    rf_key = None
    for k in history_region.historyOutputs.keys():
        ku = k.upper()
        if mode == 'circum':
            if ku == 'U2':
                u_key = k
            elif ku == 'RF2':
                rf_key = k
        elif mode == 'long':
            if ku == 'U3':
                u_key = k
            elif ku == 'RF3':
                rf_key = k

    if u_key is None or rf_key is None:
        raise RuntimeError('Could not find displacement/reaction outputs.')

    u_data = history_region.historyOutputs[u_key].data
    rf_data = history_region.historyOutputs[rf_key].data

    out = []
    for (tu, u), (tr, rf) in zip(u_data, rf_data):
        lam = 1.0 + u
        sig = abs(rf)
        out.append((lam, sig))

    fname = 'circum_job_response.txt' if mode == 'circum' else 'long_job_response.txt'

    with open(fname, 'w') as f:
        for lam, sig in out:
            f.write('%.8f %.8f\n' % (lam, sig))

    odb.close()
    print('Wrote', fname)


def _read_two_columns(path):
    xs, ys = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.replace(',', ' ').replace('\t', ' ').split()
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
    return xs, ys


def plot(response_file):
    import matplotlib.pyplot as plt

    lam, sig = _read_two_columns(response_file)

    plt.figure(figsize=(6, 5))
    plt.plot(lam, sig, 'o-', linewidth=2, markersize=5)
    plt.xlabel('Stretch (lambda)')
    plt.ylabel('Stress')
    plt.title('Simulated stress-stretch response')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('stress_strain.png', dpi=150)
    print('Wrote stress_strain.png')


def _linspace(start, stop, step):
    n = int(round((stop - start) / step))
    return [start + i * step for i in range(n + 1)]


def _interp1_linear_extrap(x_known, y_known, x_query):
    """Reimplementation of MATLAB's interp1(...,'linear','extrap')."""
    pairs = sorted(zip(x_known, y_known))
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    n = len(xs)

    out = []
    for xq in x_query:
        if xq <= xs[0]:
            x0, x1, y0, y1 = xs[0], xs[1], ys[0], ys[1]
        elif xq >= xs[-1]:
            x0, x1, y0, y1 = xs[-2], xs[-1], ys[-2], ys[-1]
        else:
            x0 = x1 = y0 = y1 = None
            for i in range(n - 1):
                if xs[i] <= xq <= xs[i + 1]:
                    x0, x1, y0, y1 = xs[i], xs[i + 1], ys[i], ys[i + 1]
                    break
        t = (xq - x0) / (x1 - x0)
        out.append(y0 + t * (y1 - y0))
    return out


def _rmse_nrmse(sigma_fit, sigma_exp):
    n = len(sigma_exp)
    sq_err = sum((f - e) ** 2 for f, e in zip(sigma_fit, sigma_exp)) / n
    rmse = sq_err ** 0.5
    nrmse = rmse / max(max(sigma_exp), 1.0)
    return rmse, nrmse


def score(response_file, exp_file, direction):
    """Full circumi/longi-style scoring + publication-style plot."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    if direction not in ('circum', 'long'):
        raise ValueError("direction must be 'circum' or 'long'")

    lambda_sim, sigma_sim = _read_two_columns(response_file)
    lambda_exp_raw, sigma_exp_raw = _read_two_columns(exp_file)

    # Resample experimental data at 0.05 stretch increments
    lambda_grid = _linspace(min(lambda_exp_raw), max(lambda_exp_raw), 0.05)
    sigma_exp_grid = _interp1_linear_extrap(lambda_exp_raw, sigma_exp_raw, lambda_grid)

    # Interpolate simulation onto the same grid
    sigma_fit_grid = _interp1_linear_extrap(lambda_sim, sigma_sim, lambda_grid)

    rmse, nrmse = _rmse_nrmse(sigma_fit_grid, sigma_exp_grid)

    tag = 'circumi' if direction == 'circum' else 'longi'
    export_name = 'passive_circum_export.txt' if direction == 'circum' else 'passive_longi_export.txt'
    metrics_name = '%s_fit_metrics.txt' % tag
    plot_name = '%s_publication_fit.png' % ('circum' if direction == 'circum' else 'long')
    muscle_label = 'circumferential smooth muscle' if direction == 'circum' else 'longitudinal smooth muscle'
    series_label = 'circumi' if direction == 'circum' else 'longi'

    with open(export_name, 'w') as f:
        for lam, se, sf in zip(lambda_grid, sigma_exp_grid, sigma_fit_grid):
            f.write('%.6f\t%.6f\t%.6f\n' % (lam, se, sf))
    print('Exported %s data to %s' % (direction, export_name))

    print('\n====================================')
    print('%s FIT METRICS (0.05 grid)' % tag.upper())
    print('RMSE  = %.10f' % rmse)
    print('NRMSE = %.10f' % nrmse)
    print('====================================')

    with open(metrics_name, 'w') as f:
        f.write('%s_RMSE %.10f\n' % (tag, rmse))
        f.write('%s_NRMSE %.10f\n' % (tag, nrmse))
    print('Wrote', metrics_name)

    # ---- publication-style figure (matches the old "Plot 2" MATLAB style) ----
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']
    mpl.rcParams['mathtext.fontset'] = 'stix'

    fig, ax = plt.subplots(figsize=(9, 7.2))

    ax.plot(lambda_grid, sigma_exp_grid, 'o',
            markersize=15,
            markerfacecolor=(0.8, 0.8, 0.8),
            markeredgecolor='k',
            markeredgewidth=3.2,
            linestyle='None',
            label='Data')

    ax.plot(lambda_grid, sigma_fit_grid, '-',
            color='k',
            linewidth=6.0,
            label='Fitted Model')

    ax.set_xlabel(r'Stretch $\lambda$', fontsize=36)
    ax.set_ylabel(r'Stress $\sigma^{\mathrm{pas}}$ (kPa)', fontsize=36)
    ax.set_title('Passive response of\n%s' % muscle_label, fontsize=36)

    ax.legend(fontsize=20, loc='upper left', frameon=True)

    ax.set_box_aspect(1)
    ax.grid(False)
    ax.set_xlim(min(lambda_grid), max(lambda_grid))

    for spine in ax.spines.values():
        spine.set_linewidth(2.4)
    ax.tick_params(width=2.4, length=8, labelsize=32)

    plt.tight_layout()
    plt.savefig(plot_name, dpi=200)
    print('Wrote', plot_name)


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if '--extract' in sys.argv:
        extract(args[0], args[1])
    elif '--score' in sys.argv:
        score(args[0], args[1], args[2])
    elif '--plot' in sys.argv:
        plot(args[0])
    else:
        print(__doc__)
        sys.exit(1)