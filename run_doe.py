# run_doe.py
#
# Driver for the full-factorial 2^5 design of experiments
# (eta_D, eta_a, eta_T, eta_mu, gamma_h).
#
# Runs the two-model workflow once per DOE combination in its own working
# directory and collects the motility metric MM into a table.
#
# Run:
#   python3 run_doe.py                  # Model 1 references + all 32 runs
#   python3 run_doe.py --dry-run        # prepare everything, launch nothing
#   python3 run_doe.py --only run05     # single case, for debugging
#   python3 run_doe.py --results-only   # re-collect MM from existing runs
#
# Per case:
#   Model 1  -> m1.odb                  (gamma_h selects the reference)
#   export   -> triad_fields.inp
#   Model 2  -> m2.odb
#   volume   -> motility_volume_cross_sections.csv
#   metric   -> MM = sigma_V over the last peristaltic cycle

from __future__ import print_function

import argparse
import os
import re
import shutil
import subprocess
import sys
import threading
import time

from concurrent.futures import ThreadPoolExecutor

import generate_doe


# ----------------------------------------------------------------------
# Layout
# ----------------------------------------------------------------------

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))

WORK_DIR = 'doe_work'
MODEL1_SUBDIR = 'model_1'
RUNS_SUBDIR = 'runs'

MODEL1_JOB = 'm1'
MODEL2_JOB = 'm2'

TRIAD_FILE = 'triad_fields.inp'
VOLUME_CSV = 'motility_volume_cross_sections.csv'

# Final DOE table, written next to run_doe.py.
SUMMARY_TXT = 'DOE_MM_results.txt'

# Files every working directory needs but never modifies.  These are
# linked rather than copied; mesh.inp alone is ~1.8 MB per case.
SHARED_LINKS_MODEL1 = ['mesh.inp', 'sub1.f', 'export_triad_fields.py']
SHARED_LINKS_MODEL2 = [
    'mesh.inp',
    'sub2.f',
    'motility_from_cross_sections.py',
    'motility_metric.py',
]

# Restart/state files Abaqus needs in the working directory.
MODEL1_IMPORT_SUFFIXES = ['.odb', '.mdl', '.prt', '.res', '.stt', '.sim',
                          '.stt.pre', '.abq', '.pac', '.sel']

# The fibrotic seed condition in m1.inp.  Both references hold it at one;
# they differ only in gamma_h.
SEED_LINE_PATTERN = re.compile(
    r'^(\s*FIBROTIC_SEED\s*,\s*11\s*,\s*11\s*,\s*)(\S+)\s*$',
    re.IGNORECASE
)

# gamma_h (Model 1's GAMMA_G, PROPS(4)) is located by the
# MECHANICAL *USER MATERIAL header on the preceding line.
GAMMA_HEADER_PATTERN = re.compile(
    r'^\s*\*USER MATERIAL\s*,\s*TYPE\s*=\s*MECHANICAL', re.IGNORECASE
)

# m2.inp's baseline *PARAMETER block (manual-run eta values).  The
# automation replaces it with an *INCLUDE of the case's own eta_params.inp.
PARAMETER_HEADER_PATTERN = re.compile(r'^\s*\*PARAMETER\s*$', re.IGNORECASE)
ETA_ASSIGNMENT_PATTERN = re.compile(r'^\s*eta_\w+\s*=\s*\S+\s*$',
                                    re.IGNORECASE)
ETA_INCLUDE_LOCAL = 'eta_params.inp'

# Emitted by Abaqus into the .sta file when an analysis finishes cleanly.
SUCCESS_MARKER = 'THE ANALYSIS HAS COMPLETED SUCCESSFULLY'

# Printed by motility_metric.py.
MM_PATTERN = re.compile(r'Motility metric,\s*MM\s*=\s*sigma_V\s*:\s*(\S+)')

FAILED = 'FAILED'
SKIPPED = 'SKIPPED'
DRY_RUN = 'DRY_RUN'


# ----------------------------------------------------------------------
# Output buffering
# ----------------------------------------------------------------------

# When cases run concurrently their log lines would interleave, so each
# worker thread collects its own lines and the whole block is printed at
# once.  With --jobs 1 there is no buffer and output streams directly.
_OUTPUT_LOCK = threading.Lock()
_THREAD_STATE = threading.local()


def log(message):
    buffer = getattr(_THREAD_STATE, 'buffer', None)
    if buffer is None:
        print(message)
        sys.stdout.flush()
    else:
        buffer.append(message)


class BufferedOutput(object):
    """Collect this thread's log lines and emit them as one block."""

    def __init__(self, header):
        self.header = header

    def __enter__(self):
        _THREAD_STATE.buffer = []
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        lines = getattr(_THREAD_STATE, 'buffer', [])
        _THREAD_STATE.buffer = None
        with _OUTPUT_LOCK:
            if self.header is not None:
                print(self.header)
            for line in lines:
                print(line)
            sys.stdout.flush()
        return False


def run_in_parallel(items, task, concurrency, header=None):
    """Apply task to each item, at most `concurrency` running at once.

    Returns a list of (item, result) in the original item order.
    """
    results = [None] * len(items)

    def worker(pair):
        index, item = pair
        title = header(item) if header else None
        if concurrency == 1:
            if title is not None:
                log(title)
            results[index] = task(item)
        else:
            with BufferedOutput(title):
                results[index] = task(item)

    pairs = list(enumerate(items))
    if concurrency == 1 or len(pairs) <= 1:
        for pair in pairs:
            worker(pair)
    else:
        with ThreadPoolExecutor(max_workers=concurrency) as pool:
            list(pool.map(worker, pairs))

    return list(zip(items, results))


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def link_or_copy(source, destination):
    if os.path.lexists(destination):
        os.remove(destination)
    try:
        os.symlink(source, destination)
    except (OSError, AttributeError, NotImplementedError):
        shutil.copy2(source, destination)


def read_text(path):
    with open(path, 'r') as handle:
        return handle.read()


def rewrite_lines(source_text, transform):
    """Apply transform to each line, preserving the file's final newline.

    Returns (text, number_of_substitutions).
    """
    output_lines = []
    substitutions = 0
    for line in source_text.splitlines():
        replacement = transform(line)
        if replacement is None:
            output_lines.append(line)
        else:
            output_lines.append(replacement)
            substitutions += 1

    text = '\n'.join(output_lines)
    if source_text.endswith('\n'):
        text += '\n'
    return text, substitutions


def write_text(path, text):
    with open(path, 'w') as handle:
        handle.write(text)


def tail(path, number_of_lines=25):
    """Return the last lines of a file, or '' if it cannot be read."""
    if not os.path.exists(path):
        return ''
    try:
        with open(path, 'r') as handle:
            lines = handle.readlines()
    except IOError:
        return ''
    return ''.join(lines[-number_of_lines:]).strip()


# ----------------------------------------------------------------------
# Deck preparation
# ----------------------------------------------------------------------

def write_model1_deck(destination, seed_value, gamma_value):
    """Copy m1.inp, setting the fibrotic seed and growth rate gamma_h."""
    source_text = read_text(os.path.join(REPO_ROOT, 'm1.inp'))

    def set_seed(line):
        match = SEED_LINE_PATTERN.match(line)
        if not match:
            return None
        return '%s%s' % (match.group(1), seed_value)

    text, substitutions = rewrite_lines(source_text, set_seed)

    if substitutions != 1:
        raise RuntimeError(
            'Expected exactly one FIBROTIC_SEED line in m1.inp, found %d. '
            'The deck has changed; run_doe.py needs updating.'
            % substitutions
        )

    # gamma_h is the 4th of the 6 comma-separated MECHANICAL constants, on
    # the line right after the TYPE=MECHANICAL header.
    pending_gamma_line = [False]

    def set_gamma(line):
        if pending_gamma_line[0]:
            pending_gamma_line[0] = False
            parts = [part.strip() for part in line.split(',')]
            if len(parts) != 6:
                return None
            parts[3] = gamma_value
            return ', '.join(parts)
        if GAMMA_HEADER_PATTERN.match(line):
            pending_gamma_line[0] = True
        return None

    text, substitutions = rewrite_lines(text, set_gamma)

    if substitutions != 1:
        raise RuntimeError(
            'Expected exactly one MECHANICAL constants line after the '
            'TYPE=MECHANICAL header in m1.inp, found %d. '
            'The deck has changed; run_doe.py needs updating.'
            % substitutions
        )

    write_text(destination, text)


def write_model2_deck(destination):
    """Copy m2.inp, replacing its baseline *PARAMETER block with an
    *INCLUDE of eta_params.inp, which prepare_case() writes alongside it.
    """
    source_text = read_text(os.path.join(REPO_ROOT, 'm2.inp'))
    lines = source_text.splitlines()

    header_indices = [index for index, line in enumerate(lines)
                      if PARAMETER_HEADER_PATTERN.match(line)]
    if len(header_indices) != 1:
        raise RuntimeError(
            'Expected exactly one *PARAMETER block in m2.inp, found %d. '
            'The deck has changed; run_doe.py needs updating.'
            % len(header_indices)
        )

    start = header_indices[0]
    end = start + 1
    while end < len(lines) and ETA_ASSIGNMENT_PATTERN.match(lines[end]):
        end += 1
    if end == start + 1:
        raise RuntimeError(
            'Found *PARAMETER in m2.inp but no eta_* assignment lines '
            'after it. The deck has changed; run_doe.py needs updating.'
        )

    new_lines = (lines[:start]
                + ['*INCLUDE, INPUT=%s' % ETA_INCLUDE_LOCAL]
                + lines[end:])
    text = '\n'.join(new_lines)
    if source_text.endswith('\n'):
        text += '\n'

    write_text(destination, text)


# ----------------------------------------------------------------------
# Command execution
# ----------------------------------------------------------------------

class Runner(object):
    """Runs the workflow commands, or reports them in dry-run."""

    def __init__(self, abaqus, python, cpus, dry_run,
                 stream_output=True):
        self.abaqus = abaqus
        self.python = python
        self.cpus = cpus
        self.dry_run = dry_run
        # With concurrent cases the Abaqus launcher output of several jobs
        # would interleave on the terminal, so it goes to a per-case file
        # instead.  Running one case at a time streams it directly.
        self.stream_output = stream_output

    def execute(self, command, cwd, capture=False, extra_env=None,
                log_path=None):
        """Return (ok, output). In dry-run nothing is launched."""
        if self.dry_run:
            log('        [dry-run] %s' % ' '.join(command))
            return True, ''

        environment = os.environ.copy()
        if extra_env:
            environment.update(extra_env)

        try:
            if capture:
                process = subprocess.Popen(
                    command, cwd=cwd, env=environment,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True
                )
                output, _ = process.communicate()
                return process.returncode == 0, output
            if log_path is not None:
                with open(os.path.join(cwd, log_path), 'w') as handle:
                    returncode = subprocess.call(
                        command, cwd=cwd, env=environment,
                        stdout=handle, stderr=subprocess.STDOUT
                    )
                return returncode == 0, ''
            returncode = subprocess.call(command, cwd=cwd, env=environment)
            return returncode == 0, ''
        except OSError as error:
            return False, 'Could not launch %s: %s' % (command[0], error)

    def abaqus_job(self, job, deck, user_subroutine, cwd, oldjob=None):
        command = [
            self.abaqus,
            'job=%s' % job,
            'input=%s' % deck,
            'user=%s' % user_subroutine,
            'cpus=%d' % self.cpus,
            'interactive',
        ]
        if oldjob:
            command.insert(3, 'oldjob=%s' % oldjob)
        # Jobs stay 'interactive'; concurrency comes from running several
        # cases at once, not from detaching individual Abaqus jobs.
        log_path = None if self.stream_output else '%s.launch.log' % job
        return self.execute(command, cwd, log_path=log_path)

    def abaqus_python(self, script, arguments, cwd):
        command = [self.abaqus, 'python', script] + arguments
        return self.execute(command, cwd, capture=True)

    def abaqus_viewer(self, script, cwd):
        command = [self.abaqus, 'viewer', 'noGUI=%s' % script]
        return self.execute(command, cwd, capture=True)

    def metric_script(self, cwd):
        return self.execute([self.python, 'motility_metric.py'], cwd,
                            capture=True)


def job_completed(directory, job):
    """True when Abaqus reported a clean finish for this job."""
    status_path = os.path.join(directory, '%s.sta' % job)
    if not os.path.exists(status_path):
        return False
    return SUCCESS_MARKER in read_text(status_path).upper()


def job_diagnostics(directory, job):
    """Return the most useful lines Abaqus left behind after a failure."""
    for suffix in ('.sta', '.msg', '.dat'):
        text = tail(os.path.join(directory, '%s%s' % (job, suffix)))
        if text:
            return '%s%s: %s' % (job, suffix, text.replace('\n', ' | '))
    return 'no Abaqus diagnostic files found'


# ----------------------------------------------------------------------
# Reference configurations (Model 1)
# ----------------------------------------------------------------------

def prepare_reference(reference_dir, seed_value, gamma_value):
    ensure_dir(reference_dir)
    for name in SHARED_LINKS_MODEL1:
        link_or_copy(os.path.join(REPO_ROOT, name),
                     os.path.join(reference_dir, name))
    write_model1_deck(os.path.join(reference_dir, 'm1.inp'),
                      seed_value, gamma_value)


def build_reference(name, reference_dir, seed_value, gamma_value, runner):
    """Run Model 1 and export the triad fields for one reference case."""
    log('  Reference "%s"  (FIBROTIC_SEED = %s, gamma_h = %s)'
        % (name, seed_value, gamma_value))
    prepare_reference(reference_dir, seed_value, gamma_value)

    triad_path = os.path.join(reference_dir, TRIAD_FILE)

    ok, _ = runner.abaqus_job(MODEL1_JOB, 'm1.inp', 'sub1.f', reference_dir)
    if not runner.dry_run:
        ok = ok and job_completed(reference_dir, MODEL1_JOB)
    if not ok:
        log('    -> Model 1: FAILED')
        log('       %s' % job_diagnostics(reference_dir, MODEL1_JOB))
        return False, 'Model 1 failed for reference %s' % name
    log('    -> Model 1: SUCCESS')

    ok, output = runner.abaqus_python(
        'export_triad_fields.py',
        ['--odb', '%s.odb' % MODEL1_JOB, '--out', TRIAD_FILE],
        reference_dir
    )
    if not runner.dry_run:
        ok = ok and os.path.exists(triad_path)
    if not ok:
        log('    -> Reference state: FAILED')
        log('       %s' % (output or 'triad export produced no output').strip())
        return False, 'Triad export failed for reference %s' % name
    log('    -> Reference state: SUCCESS')

    return True, ''


# ----------------------------------------------------------------------
# DOE cases (Model 2 and post-processing)
# ----------------------------------------------------------------------

def prepare_case(case_dir, reference_dir, run_number, total_runs, run,
                 runner):
    ensure_dir(case_dir)

    for name in SHARED_LINKS_MODEL2:
        link_or_copy(os.path.join(REPO_ROOT, name),
                     os.path.join(case_dir, name))

    # Reference configuration: predefined fields plus the restart/state
    # files that *IMPORT reads via oldjob=m1.
    link_or_copy(os.path.join(reference_dir, TRIAD_FILE),
                 os.path.join(case_dir, TRIAD_FILE))
    for suffix in MODEL1_IMPORT_SUFFIXES:
        source = os.path.join(reference_dir, '%s%s' % (MODEL1_JOB, suffix))
        if os.path.exists(source):
            link_or_copy(source,
                         os.path.join(case_dir,
                                      '%s%s' % (MODEL1_JOB, suffix)))

    eta_text = generate_doe.render_include(run_number, total_runs, run)
    write_text(os.path.join(case_dir, ETA_INCLUDE_LOCAL), eta_text)
    write_model2_deck(os.path.join(case_dir, 'm2.inp'))


def extract_mm(case_dir, runner):
    """Run motility_metric.py and read MM from its output."""
    ok, output = runner.metric_script(case_dir)
    if runner.dry_run:
        return None, ''
    if not ok:
        return None, (output or 'motility_metric.py produced no output')[-800:]

    match = MM_PATTERN.search(output)
    if not match:
        return None, 'MM line not found in motility_metric.py output'
    try:
        return float(match.group(1)), ''
    except ValueError:
        return None, 'Could not parse MM value %r' % match.group(1)


def run_case(case, runner):
    """Execute one DOE case end to end.  Returns the updated case dict."""
    case_dir = case['dir']
    prepare_case(case_dir, case['reference_dir'], case['run_number'],
                case['total_runs'], case['run'], runner)

    if runner.dry_run:
        case['status'] = DRY_RUN
        log('    -> prepared %s' % os.path.relpath(case_dir, REPO_ROOT))
        return case

    ok, _ = runner.abaqus_job(MODEL2_JOB, 'm2.inp', 'sub2.f', case_dir,
                              oldjob=MODEL1_JOB)
    ok = ok and job_completed(case_dir, MODEL2_JOB)
    if not ok:
        log('    -> Model 2: FAILED')
        log('       %s' % job_diagnostics(case_dir, MODEL2_JOB))
        case['status'] = FAILED
        case['error'] = job_diagnostics(case_dir, MODEL2_JOB)
        return case
    log('    -> Model 2: SUCCESS')

    ok, output = runner.abaqus_viewer('motility_from_cross_sections.py',
                                      case_dir)
    ok = ok and os.path.exists(os.path.join(case_dir, VOLUME_CSV))
    if not ok:
        log('    -> Lumen volume: FAILED')
        log('       %s' % (output or 'no output').strip()[-500:])
        case['status'] = FAILED
        case['error'] = 'motility_from_cross_sections.py failed'
        return case
    log('    -> Lumen volume: SUCCESS')

    mm_value, error = extract_mm(case_dir, runner)
    if mm_value is None:
        log('    -> MM: FAILED')
        log('       %s' % error)
        case['status'] = FAILED
        case['error'] = error
        return case

    case['MM'] = mm_value
    case['status'] = 'SUCCESS'
    log('    -> MM: %.6f' % mm_value)
    return case


# ----------------------------------------------------------------------
# Case construction
# ----------------------------------------------------------------------

def build_cases(work_dir):
    """Return the 32 DOE cases in manuscript order."""
    model_1_dir = os.path.join(work_dir, MODEL1_SUBDIR)
    runs_dir = os.path.join(work_dir, RUNS_SUBDIR)

    # Two Model 1 cases, 'baseline' leaves (theta_h = 1)
    # the geometry undeformed, 'fibrotic' grows (theta_h ~2).
    baseline_ref = os.path.join(model_1_dir, 'baseline')
    fibrotic_ref = os.path.join(model_1_dir, 'fibrotic')

    runs = generate_doe.doe_table()
    total_runs = len(runs)
    cases = []

    for index, run in enumerate(runs):
        run_number = index + 1
        name = 'run%02d' % run_number
        levels = dict((factor, level) for factor, level, _ in run)
        reference_dir = (fibrotic_ref if levels['gamma_h'] == 'M'
                         else baseline_ref)
        cases.append({
            'name': name,
            'label': ', '.join('%s=%s' % (factor, level)
                               for factor, level, _ in run),
            'run_number': run_number,
            'total_runs': total_runs,
            'run': run,
            'levels': levels,
            'values': dict((factor, value) for factor, _, value in run),
            'dir': os.path.join(runs_dir, name),
            'reference_dir': reference_dir,
            'status': 'PENDING',
            'MM': None,
            'error': '',
        })

    return cases, baseline_ref, fibrotic_ref


# ----------------------------------------------------------------------
# Result collection
# ----------------------------------------------------------------------

def format_mm(case):
    """Render MM, never as a number for a case that did not succeed."""
    if case['MM'] is not None:
        return '%.6f' % case['MM']
    if case['status'] in (DRY_RUN, SKIPPED):
        return case['status']
    return FAILED


def write_summary(cases):
    """Write the DOE table: factor levels, MM, and dMM relative to run01."""
    factor_names = generate_doe.factor_names()

    # Run 01 (all factors at baseline) is the reference for dMM (Eq. 46).
    run01 = [case for case in cases if case['name'] == 'run01']
    mm0 = run01[0]['MM'] if run01 else None

    lines = []
    lines.append('Full-factorial 2^%d DOE results' % len(factor_names))
    lines.append('')

    header = '%4s  %-3s %-3s %-3s %-4s %-4s  %12s  %9s' % (
        'Run', 'e_D', 'e_a', 'e_T', 'e_mu', 'g_h', 'MM', 'dMM(%)'
    )
    lines.append(header)
    lines.append('-' * len(header))

    for case in cases:
        if case['MM'] is None or mm0 in (None, 0.0):
            delta_text = '%9s' % '-'
        else:
            delta_text = '%+9.1f' % ((case['MM'] - mm0) / mm0 * 100.0)
        lines.append('%4d  %-3s %-3s %-3s %-4s %-4s  %12s  %9s' % (
            case['run_number'],
            case['levels']['eta_D'], case['levels']['eta_a'],
            case['levels']['eta_T'], case['levels']['eta_mu'],
            case['levels']['gamma_h'],
            format_mm(case), delta_text
        ))

    lines.append('')
    failures = [case for case in cases if case['status'] == FAILED]
    skipped = [case['name'] for case in cases if case['status'] == SKIPPED]
    if failures:
        lines.append('FAILED cases:')
        for case in failures:
            lines.append('  %s: %s' % (case['name'], case['error']))
    else:
        lines.append('No failed cases.')
    if skipped:
        lines.append('SKIPPED cases (never started): %s' % ', '.join(skipped))

    text = '\n'.join(lines) + '\n'
    summary_path = os.path.join(REPO_ROOT, SUMMARY_TXT)
    write_text(summary_path, text)
    return text, summary_path


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Run the full-factorial 2^5 DOE using the two-model '
                    'Abaqus workflow.'
    )
    parser.add_argument('--work-dir', default=WORK_DIR,
                        help='Directory holding the Model 1 references and '
                             'the per-case run directories')
    parser.add_argument('--abaqus', default='abaqus',
                        help='Abaqus executable')
    parser.add_argument('--python', default=sys.executable,
                        help='Python used for motility_metric.py '
                             '(needs numpy)')
    parser.add_argument('--cpus', type=int, default=4,
                        help='CPUs per Abaqus job')
    parser.add_argument('--jobs', type=int, default=1,
                        help='DOE cases to run concurrently (default 1). '
                             'Peak demand is --jobs x --cpus cores, and the '
                             'Abaqus licence draw scales the same way.')
    parser.add_argument('--only', default=None,
                        help='Run a single case, e.g. run05')
    parser.add_argument('--dry-run', action='store_true',
                        help='Prepare directories and decks, launch nothing')
    parser.add_argument('--skip-model-1', action='store_true',
                        help='Reuse existing Model 1 reference runs')
    parser.add_argument('--results-only', action='store_true',
                        help='Re-extract MM from existing run directories')
    parser.add_argument('--stop-on-error', action='store_true',
                        help='Stop at the first failed case')
    return parser.parse_args()


def main():
    args = parse_args()

    work_dir = args.work_dir
    if not os.path.isabs(work_dir):
        work_dir = os.path.join(REPO_ROOT, work_dir)

    concurrency = max(1, args.jobs)
    runner = Runner(args.abaqus, args.python, args.cpus, args.dry_run,
                    stream_output=(concurrency == 1))
    cases, baseline_ref, fibrotic_ref = build_cases(work_dir)

    if args.only:
        selected = [case for case in cases if case['name'] == args.only]
        if not selected:
            log('Unknown case %r. Available: %s'
                % (args.only, ', '.join(case['name'] for case in cases)))
            return 2
        cases = selected

    log('')
    log('Full-factorial 2^%d DOE  (%d case%s)'
        % (len(generate_doe.factor_names()), len(cases),
           '' if len(cases) == 1 else 's'))
    log('Working directory: %s' % os.path.relpath(work_dir, REPO_ROOT))
    if concurrency > 1:
        log('Concurrency: %d cases at a time, %d CPUs each '
            '(peak %d CPUs; Abaqus licence draw scales with this).'
            % (concurrency, args.cpus, concurrency * args.cpus))
        log('Abaqus launcher output goes to <job>.launch.log in each '
            'case directory.')
    if args.dry_run:
        log('DRY RUN: decks and directories are prepared, '
            'no Abaqus job is launched.')
    log('')

    started = time.time()

    # ---- Reference configurations -----------------------------------
    if not args.results_only:
        needed = set(case['reference_dir'] for case in cases)
        gamma_b = generate_doe.GAMMA_H_LEVELS['B']
        gamma_m = generate_doe.GAMMA_H_LEVELS['M']
        reference_specs = [
            ('baseline (gamma_h=B, ungrown)', baseline_ref, '1.0',
             repr(gamma_b)),
            ('fibrotic (gamma_h=M, grown)', fibrotic_ref, '1.0',
             repr(gamma_m)),
        ]
        log('Building reference configurations (Model 1)')
        pending = []
        for name, reference_dir, seed_value, gamma_value in reference_specs:
            if reference_dir not in needed:
                continue
            if args.skip_model_1 and os.path.exists(
                    os.path.join(reference_dir, TRIAD_FILE)):
                log('  Reference "%s": reusing existing run' % name)
                continue
            pending.append((name, reference_dir, seed_value, gamma_value))

        # Every Model 2 case depends on a reference, so all references
        # must finish before any case starts.  They are independent of
        # each other and may overlap.
        outcomes = run_in_parallel(
            pending,
            lambda spec: build_reference(spec[0], spec[1], spec[2],
                                         spec[3], runner),
            concurrency
        )
        for _, (ok, error) in outcomes:
            if not ok:
                log('')
                log('ABORTING: %s' % error)
                return 1
        log('')

    # ---- Cases -------------------------------------------------------
    total = len(cases)
    case_positions = dict((case['name'], index + 1)
                          for index, case in enumerate(cases))
    # With --stop-on-error and concurrent cases, already-running cases are
    # allowed to finish; only cases not yet started are abandoned.
    abandon = threading.Event()

    def case_header(case):
        return 'Run %d/%d: %s  [%s]' % (
            case_positions[case['name']], total, case['label'], case['name']
        )

    def process_case(case):
        if abandon.is_set():
            case['status'] = SKIPPED
            case['error'] = 'not started (earlier case failed)'
            log('    -> SKIPPED (an earlier case failed)')
            return case

        if args.results_only:
            if not os.path.isdir(case['dir']):
                log('    -> directory missing, skipping')
                case['status'] = FAILED
                case['error'] = 'run directory not found'
                return case
            mm_value, error = extract_mm(case['dir'], runner)
            if mm_value is None:
                log('    -> MM: FAILED')
                log('       %s' % error)
                case['status'] = FAILED
                case['error'] = error
            else:
                case['MM'] = mm_value
                case['status'] = 'SUCCESS'
                log('    -> MM: %.6f' % mm_value)
        else:
            run_case(case, runner)

        if case['status'] == FAILED and args.stop_on_error:
            log('    -> --stop-on-error: no further cases will start.')
            abandon.set()
        return case

    run_in_parallel(cases, process_case, concurrency, header=case_header)

    # ---- Results -----------------------------------------------------
    log('')
    summary, summary_path = write_summary(cases)

    log(summary)
    log('Results  : %s' % os.path.relpath(summary_path, REPO_ROOT))
    log('Elapsed  : %.1f s' % (time.time() - started))

    incomplete = [case for case in cases
                  if case['status'] in (FAILED, SKIPPED)]
    return 1 if incomplete else 0


if __name__ == '__main__':
    sys.exit(main())
