# motility_metric.py
#
# Computes the manuscript motility metric
#     MM = sigma_V / V_bar
# over the last complete contraction cycle of V(t).

import csv
import numpy as np
import matplotlib.pyplot as plt

CSV_IN = 'motility_volume_cross_sections.csv'
TXT_OUT = 'cycle_volume.txt'

MIN_PERIOD_FACTOR = 0.6
MAX_PERIOD_FACTOR = 1.4
END_WINDOW_PERIODS = 1.5

# The final cycle is searched near the end of the simulation over this
# many dominant periods, then paired with the preceding compatible peak.


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def local_maxima(values):
    """Return indices of local maxima in a raw signal."""
    maxima_indices = []
    index = 1
    while index < len(values) - 1:
        is_local_maximum = (
            values[index] >= values[index - 1]
            and values[index] >= values[index + 1]
            and (
                values[index] > values[index - 1]
                or values[index] > values[index + 1]
            )
        )
        if is_local_maximum:
            maxima_indices.append(index)
        index += 1
    return maxima_indices


# ----------------------------------------------------------------------
# Read V(t)
# ----------------------------------------------------------------------

time_values = []
volume_values = []

with open(CSV_IN, 'r') as csv_file:
    reader = csv.DictReader(csv_file)
    for row in reader:
        time_values.append(float(row['time']))
        volume_values.append(float(row['volume']))

if len(time_values) < 10:
    raise RuntimeError('Not enough data points in %s' % CSV_IN)

# ----------------------------------------------------------------------
# Full-signal summary metrics
# ----------------------------------------------------------------------

mean_volume_full = np.mean(volume_values)
std_volume_full = np.std(volume_values, ddof=0)
volume_min_full = min(volume_values)
volume_max_full = max(volume_values)
peak_to_peak_full = volume_max_full - volume_min_full

# ----------------------------------------------------------------------
# Dominant contraction period from FFT
# ----------------------------------------------------------------------

delta_t = time_values[1] - time_values[0]
volume_signal = np.array(volume_values, dtype=float)
volume_signal = volume_signal - np.mean(volume_signal)

frequencies = np.fft.rfftfreq(len(volume_signal), d=delta_t)
fft_amplitudes = np.abs(np.fft.rfft(volume_signal))


# Ignore the zero-frequency component and identify the dominant nonzero
# frequency of the volume signal.
dominant_frequency_index = np.argmax(fft_amplitudes[1:]) + 1
dominant_frequency = frequencies[dominant_frequency_index]
dominant_period = 1.0 / dominant_frequency

# ----------------------------------------------------------------------
# Last complete cycle by peak detection
# ----------------------------------------------------------------------

peak_indices = local_maxima(volume_values)

if len(peak_indices) < 2:
    raise RuntimeError('Could not find enough peaks')

end_time_threshold = time_values[-1] - END_WINDOW_PERIODS * dominant_period

candidate_end_peaks = [
    peak_index for peak_index in peak_indices
    if time_values[peak_index] >= end_time_threshold
]
if not candidate_end_peaks:
    candidate_end_peaks = peak_indices


# Use the largest nearby peak as the end of the final representative cycle.
cycle_end_index = max(
    candidate_end_peaks,
    key=lambda peak_index: volume_values[peak_index]
)

best_start_index = None
best_period_error = None

for peak_index in peak_indices:
    if peak_index >= cycle_end_index:
        continue

    cycle_duration_candidate = (
        time_values[cycle_end_index] - time_values[peak_index]
    )

    within_period_window = (
        MIN_PERIOD_FACTOR * dominant_period
        <= cycle_duration_candidate
        <= MAX_PERIOD_FACTOR * dominant_period
    )

    if within_period_window:
        period_error = abs(cycle_duration_candidate - dominant_period)

        if (
            best_period_error is None
            or period_error < best_period_error
            or (
                abs(period_error - best_period_error) < 1e-12
                and volume_values[peak_index]
                > volume_values[best_start_index]
            )
        ):
            best_period_error = period_error
            best_start_index = peak_index


# If no peak lies inside the preferred period window, fall back to the
# preceding peak whose duration is closest to the dominant period.
if best_start_index is None:
    for peak_index in peak_indices:
        if peak_index >= cycle_end_index:
            continue

        cycle_duration_candidate = (
            time_values[cycle_end_index] - time_values[peak_index]
        )
        period_error = abs(cycle_duration_candidate - dominant_period)

        if best_period_error is None or period_error < best_period_error:
            best_period_error = period_error
            best_start_index = peak_index

cycle_start_index = best_start_index

cycle_time_values = time_values[cycle_start_index:cycle_end_index + 1]
cycle_volume_values = volume_values[cycle_start_index:cycle_end_index + 1]

# ----------------------------------------------------------------------
# Cycle metrics and manuscript motility metric
# ----------------------------------------------------------------------

mean_volume_cycle = np.mean(cycle_volume_values)
std_volume_cycle = np.std(cycle_volume_values, ddof=0)
volume_min_cycle = min(cycle_volume_values)
volume_max_cycle = max(cycle_volume_values)
peak_to_peak_cycle = volume_max_cycle - volume_min_cycle


# Manuscript motility metric: coefficient of variation of V(t) over
# the selected cycle.
MM = std_volume_cycle / mean_volume_cycle if mean_volume_cycle != 0 else 0.0
peak_to_peak_over_mean = (
    peak_to_peak_cycle / mean_volume_cycle
    if mean_volume_cycle != 0
    else 0.0
)

# ----------------------------------------------------------------------
# Print
# ----------------------------------------------------------------------

print('\n================ METRICS ================')
print('Cycle start time           :', time_values[cycle_start_index])
print('Cycle end time             :', time_values[cycle_end_index])
print('Cycle duration             :',
      time_values[cycle_end_index] - time_values[cycle_start_index])
print('Mean volume (cycle)        :', mean_volume_cycle)
print('STD volume (cycle)         :', std_volume_cycle)
print('STD / mean  (MM)           :', MM)
print('Peak-to-peak               :', peak_to_peak_cycle)
print('Peak-to-peak / mean        :', peak_to_peak_over_mean)
print('========================================\n')

# ----------------------------------------------------------------------
# Export cycle data
# ----------------------------------------------------------------------

with open(TXT_OUT, 'w') as text_file:
    text_file.write('time\tvolume\n')
    for time_value, volume_value in zip(cycle_time_values,
                                        cycle_volume_values):
        text_file.write('%.6f\t%.6f\n' % (time_value, volume_value))

print('Cycle data exported to %s  (%d points)\n'
      % (TXT_OUT, len(cycle_time_values)))

# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------

plt.figure(figsize=(9, 4.5))
plt.plot(time_values, volume_values, label='Volume')
plt.axvspan(time_values[cycle_start_index],
            time_values[cycle_end_index],
            alpha=0.3, label='Last cycle')
plt.scatter([time_values[cycle_start_index], time_values[cycle_end_index]],
            [volume_values[cycle_start_index],
             volume_values[cycle_end_index]],
            color='red', label='Peaks')

plt.xlabel('Time')
plt.ylabel('Volume')
plt.title('Lumen Volume vs Time')
plt.legend()
plt.grid(True)

plt.show()
