# report_motility_metrics.py

import csv
import numpy as np

CSV_IN = 'motility_volume_cross_sections.csv'

MIN_PERIOD_FACTOR = 0.6
MAX_PERIOD_FACTOR = 1.4
END_WINDOW_PERIODS = 1.5

# Number of equally spaced points used to calculate MM
N_METRIC_POINTS = 1000


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def local_minima(vals):
    idx = []
    i = 1

    while i < len(vals) - 1:
        if vals[i] <= vals[i - 1] and vals[i] <= vals[i + 1] and (
            vals[i] < vals[i - 1] or vals[i] < vals[i + 1]
        ):
            idx.append(i)

        i += 1

    return idx


# ----------------------------------------------------------------------
# Read CSV
# ----------------------------------------------------------------------

time = []
volume = []

with open(CSV_IN, 'r') as f:
    reader = csv.DictReader(f)

    for row in reader:
        time.append(float(row['time']))
        volume.append(float(row['volume']))

if len(time) < 10:
    raise RuntimeError('Not enough data points in %s' % CSV_IN)


# ----------------------------------------------------------------------
# Convert to arrays
# ----------------------------------------------------------------------

time = np.array(time, dtype=float)
volume = np.array(volume, dtype=float)


# ----------------------------------------------------------------------
# Check that time is increasing
# ----------------------------------------------------------------------

if np.any(np.diff(time) <= 0):
    raise RuntimeError('Time values must be strictly increasing')


# ----------------------------------------------------------------------
# Dominant period
#
# FFT requires equally spaced data, so interpolate the full signal first.
# ----------------------------------------------------------------------

dt_raw = np.diff(time)
dt_fft = np.mean(dt_raw)

time_fft = np.arange(
    time[0],
    time[-1] + 0.5 * dt_fft,
    dt_fft
)

volume_fft = np.interp(
    time_fft,
    time,
    volume
)

x = volume_fft - np.mean(volume_fft)

freq = np.fft.rfftfreq(
    len(x),
    d=dt_fft
)

fft_vals = np.abs(
    np.fft.rfft(x)
)

if len(fft_vals) < 2:
    raise RuntimeError('Not enough data for FFT')

# Ignore zero-frequency component
best_k = np.argmax(fft_vals[1:]) + 1
f_dom = freq[best_k]

if f_dom <= 0:
    raise RuntimeError('Could not determine a valid dominant frequency')

T_dom = 1.0 / f_dom


# ----------------------------------------------------------------------
# Valley detection on raw signal
# ----------------------------------------------------------------------

valleys = local_minima(volume)

if len(valleys) < 2:
    raise RuntimeError('Could not find enough valleys')


# ----------------------------------------------------------------------
# Choose ending valley near end of simulation
# ----------------------------------------------------------------------

t_end_threshold = time[-1] - END_WINDOW_PERIODS * T_dom

candidate_end = [
    i for i in valleys
    if time[i] >= t_end_threshold
]

if not candidate_end:
    candidate_end = valleys

i_end = min(
    candidate_end,
    key=lambda i: volume[i]
)


# ----------------------------------------------------------------------
# Find starting valley approximately one dominant period earlier
# ----------------------------------------------------------------------

best_start = None
best_err = None

for idx in valleys:

    if idx >= i_end:
        continue

    dt_cycle = time[i_end] - time[idx]

    if MIN_PERIOD_FACTOR * T_dom <= dt_cycle <= MAX_PERIOD_FACTOR * T_dom:

        err = abs(dt_cycle - T_dom)

        if (
            best_err is None
            or err < best_err
            or (
                abs(err - best_err) < 1e-12
                and volume[idx] < volume[best_start]
            )
        ):
            best_err = err
            best_start = idx


# ----------------------------------------------------------------------
# Fallback if no start valley lies within allowed period range
# ----------------------------------------------------------------------

if best_start is None:

    for idx in valleys:

        if idx >= i_end:
            continue

        dt_cycle = time[i_end] - time[idx]
        err = abs(dt_cycle - T_dom)

        if best_err is None or err < best_err:
            best_err = err
            best_start = idx


if best_start is None:
    raise RuntimeError('Could not identify starting valley')


i_start = best_start


# ----------------------------------------------------------------------
# Extract raw selected cycle
# ----------------------------------------------------------------------

cycle_time_raw = time[i_start:i_end + 1]
cycle_volume_raw = volume[i_start:i_end + 1]

if len(cycle_volume_raw) < 2:
    raise RuntimeError('Selected cycle contains too few points')


# ----------------------------------------------------------------------
# Interpolate selected cycle onto equally spaced grid
# ----------------------------------------------------------------------

cycle_time = np.linspace(
    cycle_time_raw[0],
    cycle_time_raw[-1],
    N_METRIC_POINTS
)

cycle_volume = np.interp(
    cycle_time,
    cycle_time_raw,
    cycle_volume_raw
)


# ----------------------------------------------------------------------
# Motility metric
#
# MM = standard deviation of lumen volume over the selected cycle
# ----------------------------------------------------------------------

MM = np.std(cycle_volume)

cycle_duration = cycle_time_raw[-1] - cycle_time_raw[0]


# ----------------------------------------------------------------------
# Print
# ----------------------------------------------------------------------


print('Motility metric, MM = sigma_V :', MM)
print('Cycle duration                :', cycle_duration)
