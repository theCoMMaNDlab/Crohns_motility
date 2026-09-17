# Crohn's Motility — Computational Study of Fibrotic Strictures on Intestinal Motility

A two-model Abaqus finite element framework simulating the effect of fibrotic strictures on intestinal electromechanics and peristaltic motility in Crohn's disease.

---

## Overview

The framework consists of two sequential Abaqus jobs:

| Model                             | Subroutine    | Purpose                                                                                                                        |
| --------------------------------- | ------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| **Model 1** (`m1.inp` + `sub1.f`) | UMAT / UMATHT | Fibrosis field evolution (logistic reaction–diffusion) and radial wall growth via multiplicative decomposition **F = Fᵉ · Fʰ** |
| **Model 2** (`m2.inp` + `sub2.f`) | UMAT / UMATHT | FitzHugh–Nagumo electrophysiology, SMC active tension, and mechanoelectrical feedback, run on the grown geometry from Model 1  |

The temperature DOF is used as a surrogate scalar field in both models to solve the associated transport equations within the Abaqus coupled temperature–displacement framework.

The end-to-end workflow has three stages:

1. **Passive constitutive calibration** (`passive-calibration/`) — fit the passive fiber stiffnesses used in Models 1 and 2 against uniaxial stress–stretch data.
2. **Abaqus simulation** (`simulation-files/`) — run Model 1 followed by Model 2 for a given disease-parameter combination, manually or automatically across the full DOE.
3. **DOE / postprocessing** (`postprocessing-files/`) — perform regression analysis over the DOE motility-metric results.

---

## Repository Structure

```text
Crohns_motility/
├── passive-calibration/
│   ├── tune_fiber_stiffness.py          # Sequential fit of k_c, k_l against uniaxial data
│   ├── circum_stretch_stress.txt        # Circumferential stress–stretch data
│   └── long_stretch_stress.txt          # Longitudinal stress–stretch data
├── simulation-files/
│   ├── mesh.inp                         # Quarter-cylinder mesh
│   ├── m1.inp                           # Model 1 input: fibrosis growth job
│   ├── m2.inp                           # Model 2 input: electromechanical job
│   ├── sub1.f                           # UMAT/UMATHT for Model 1 (growth)
│   ├── sub2.f                           # UMAT/UMATHT for Model 2 (electromechanics)
│   ├── export_triad_fields.py           # Exports fibrosis field and fiber triad from m1.odb
│   ├── motility_from_cross_sections.py  # Computes lumen volume from m2.odb
│   ├── motility_metric.py               # Computes the motility metric
│   ├── generate_doe.py                  # Defines the DOE parameters and run table
│   └── run_doe.py                       # DOE driver: runs the workflow for every combination
├── postprocessing-files/
│   ├── DOE_MM_results.txt               # Generated: DOE summary table (run_doe.py output)
│   └── regression_analysis.py           # Regression and Figure 6-style plots over DOE_MM_results.txt
├── requirements.txt                     # Python dependencies for the standalone scripts
└── LICENSE                              # GNU General Public License v3.0
```

`triad_fields.inp` and `motility_volume_cross_sections.csv` are generated at runtime in the working directory (see Workflow below).

---

## Prerequisites

* Abaqus 2024 or later (also used to run `export_triad_fields.py` and `motility_from_cross_sections.py` via `abaqus python` / `abaqus viewer`)
* Python 3.8+ with `numpy`, `scipy`, and `matplotlib`

```bash
pip install -r requirements.txt
```

`export_triad_fields.py` and `motility_from_cross_sections.py` import `odbAccess` / `abaqusConstants`, which ship with Abaqus and are not available on PyPI. Run these two scripts with the Abaqus interpreter (`abaqus python` / `abaqus viewer`) as shown below.

---

## 1. Passive constitutive calibration

`passive-calibration/tune_fiber_stiffness.py` sequentially fits the circumferential and longitudinal fiber stiffnesses `k_c` and `k_l` against the uniaxial stress–stretch data in that folder. Run it from inside `passive-calibration/`:

```bash
python3 tune_fiber_stiffness.py
```

The script writes the fitted parameters and fit-quality metrics/plots into the same folder. These outputs are regenerated on each run and are not version-controlled.

---

## 2. Abaqus simulation

Run the steps below from inside `simulation-files/`.

### Step 1 — Run Model 1 (fibrosis growth)

```bash
abaqus job=m1 input=m1.inp user=sub1.f cpus=4 interactive
```

Outputs: `m1.odb` containing the grown geometry and the fibrosis field stored in `NT11` (the temperature DOF).

---

### Step 2 — Export triad fields

```bash
abaqus python export_triad_fields.py --odb m1.odb --out triad_fields.inp
```

This reads the last frame of `m1.odb` and writes `triad_fields.inp` containing 10 predefined field variables:

| Variable | Field                              |
| -------- | ---------------------------------- |
| 1        | Fibrosis scalar *m*                |
| 2–4      | Radial unit vector **er**          |
| 5–7      | Circumferential unit vector **eθ** |
| 8–10     | Axial unit vector **ez**          |

### Step 3 — Run Model 2 (electromechanics)

```bash
abaqus job=m2 oldjob=m1 input=m2.inp user=sub2.f cpus=4 interactive
```

Model 2 imports the grown geometry from Model 1 and reads `triad_fields.inp` as predefined fields. The output is `m2.odb`.

The `*PARAMETER` block at the top of `m2.inp` sets `eta_D`, `eta_a`, `eta_T`, and `eta_mu` for the manual run (baseline = 0.0). Edit these values directly to run a different disease level. For DOE runs, `run_doe.py` writes the corresponding parameter values for each case.

---

### Step 4 — Compute lumen volume

```bash
abaqus viewer noGUI=motility_from_cross_sections.py
```

This reads `m2.odb` and `mesh.inp`, integrates lumen cross-sectional areas along the tube axis frame by frame, and writes:

`motility_volume_cross_sections.csv`

---

### Step 5 — Compute motility metric

```bash
python3 motility_metric.py
```

This reads `motility_volume_cross_sections.csv` and computes:

* Dominant contraction period (FFT)
* Last complete contraction cycle (valley detection)
* **Motility metric** MM = STD(volume) over the last cycle

The script prints the motility metric and cycle duration to the console.

---

### Quick reference

The five manual steps in order, run from `simulation-files/`:

```bash
abaqus job=m1 input=m1.inp user=sub1.f cpus=4 interactive
abaqus python export_triad_fields.py --odb m1.odb --out triad_fields.inp
abaqus job=m2 oldjob=m1 input=m2.inp user=sub2.f cpus=4 interactive
abaqus viewer noGUI=motility_from_cross_sections.py
python3 motility_metric.py
```

Or run all 32 DOE cases end to end from `simulation-files/`:

```bash
python3 run_doe.py --cpus 8
```

**Note:** The DOE runs all 32 cases sequentially by default. Use `--jobs` to run cases concurrently; total CPU usage scales approximately with `--jobs × --cpus`.

---

## 3. DOE / postprocessing

Use `run_doe.py` to orchestrate the five manual steps above for all 32 DOE combinations (`eta_D`, `eta_a`, `eta_T`, `eta_mu`, `gamma_h`).

It runs the same commands for each of the 32 cases and writes the final table to `postprocessing-files/DOE_MM_results.txt`.

Model 1 is run twice, once for the `baseline` and once for the `fibrotic` (deformed) reference state. Each reference state is shared by the 16 DOE cases that use it, and Model 2 runs once per case in `simulation-files/doe_work/runs`.

```bash
python3 run_doe.py --cpus 8
```

Each case's `*PARAMETER` values come from the factor table in `generate_doe.py` and are written directly to that case's `eta_params.inp`.

Once `postprocessing-files/DOE_MM_results.txt` has been generated, run the regression analysis and reproduce the Figure 6-style plots from `postprocessing-files/`:

```bash
python3 regression_analysis.py
```

---

## Notes

* Model 2 uses a quarter-cylinder geometry with symmetry boundary conditions on `FACEX` and `FACEY` and axial fixity on `PROXIMAL_FACE` and `DISTAL_FACE`.
* The temperature DOF carries the normalized ICC potential *φ* ∈ [0, 1] in Model 2.
* `INNER_SURFACE` must be defined in `mesh.inp` for the volume computation to work. It represents the lumen-facing nodes used by the volume calculation.
* All scripts assume that the working directory contains the relevant `.odb`, `.inp`, and `.csv` files.
* Run `python3 run_doe.py --help` for available options, including `--jobs` for concurrent cases, `--dry-run`, `--only`, and `--results-only`.

---

## License

Released under the GNU General Public License v3.0 — see [LICENSE](LICENSE).
