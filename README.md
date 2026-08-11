# Crohn's Motility — Computational Study of Fibrotic Strictures on Intestinal Motility

A two-model Abaqus finite element framework simulating the effect of fibrotic strictures on intestinal electromechanics and peristaltic motility in Crohn's disease.

---

## Overview

The framework consists of two sequential Abaqus jobs:

| Model | Subroutine | Purpose |
|---|---|---|
| **Model 1** (`m1.inp` + `sub1.f`) | UMAT / UMATHT | Fibrosis field evolution (logistic reaction–diffusion) and radial wall growth via multiplicative decomposition **F = Fᵉ · Fʰ** |
| **Model 2** (`m2.inp` + `sub2.f`) | UMAT / UMATHT | FitzHugh–Nagumo electrophysiology, SMC active tension, mechanoelectrical feedback — run on the grown geometry from Model 1 |

The temperature DOF is used as a surrogate scalar field in both models to solve the associated transport equations within the Abaqus coupled temperature–displacement framework.

---

## Repository Structure

```
├── mesh.inp                      # Quarter-cylinder mesh 
├── m1.inp                        # Model 1 input: fibrosis growth job
├── m2.inp                        # Model 2 input: electromechanical job
├── sub1.f                        # UMAT/UMATHT for Model 1 (growth)
├── sub2.f                        # UMAT/UMATHT for Model 2 (electromechanics)
├── export_triad_fields.py        # Exports fibrosis field + fibre triad from m1.odb
├── triad_fields.inp              # Generated: predefined field include for Model 2
├── motility_from_cross_sections.py  # Computes lumen volume from m2.odb
├── motility_metric.py            # Computes motility metric (STD/mean) and plots
├── motility_volume_cross_sections.csv  # Generated: volume vs time output
└── commands.txt                  # Quick reference for run commands
```

---

## Prerequisites

- Abaqus 2024 or above (also used to run `export_triad_fields.py` and `motility_from_cross_sections.py` via `abaqus python` / `abaqus viewer`)
- Python 3.8+ with `numpy` and `matplotlib`, for the standalone post-processing script `motility_metric.py`
---

## Workflow


### Step 1 — Run Model 1 (fibrosis growth)

```bash
abaqus job=m1 input=m1.inp user=sub1.f cpus=4 interactive
```

Outputs: `m1.odb` containing the grown geometry and fibrosis field `NT11`.

---

### Step 2 — Export triad fields

```bash
abaqus python export_triad_fields.py --odb m1.odb --out triad_fields.inp
```

This reads the last frame of `m1.odb` and writes `triad_fields.inp` containing 10 predefined field variables:

| Variable | Field |
|---|---|
| 1 | Fibrosis scalar *m* |
| 2–4 | Radial unit vector **eᵣ** |
| 5–7 | Circumferential unit vector **e**θ |
| 8–10 | Axial unit vector **e**_z |

### Step 3 — Run Model 2 (electromechanics)

```bash
abaqus job=m2 oldjob=m1 input=m2.inp user=sub2.f cpus=4 interactive
```

Model 2 imports the grown geometry from Model 1 and reads `triad_fields.inp` as predefined fields. Outputs: `m2.odb`.

---

### Step 4 — Compute lumen volume

```bash
abaqus viewer noGUI=motility_from_cross_sections.py
```

Reads `m2.odb` and `mesh.inp`, integrates lumen cross-sectional areas along the tube axis frame-by-frame, and writes:

```
motility_volume_cross_sections.csv   (columns: time, volume, v_rel)
```

---

### Step 5 — Compute motility metric

```bash
python motility_metric.py
```

Reads `motility_volume_cross_sections.csv` and computes:

- Dominant contraction period (FFT)
- Last complete contraction cycle (peak detection)
- **Motility metric** = STD(volume) / mean(volume) over the last cycle

Prints metrics to console and plots lumen volume vs time with the analysis window highlighted.

---


## Notes

- Model 2 uses a quarter-cylinder geometry (symmetry boundary conditions on `FACEX` and `FACEY`; axial fixity on `PROXIMAL_FACE` and `DISTAL_FACE`).
- The temperature DOF carries the normalised ICC potential *φ* ∈ [0, 1] in Model 2.
- `INNER_SURFACE` must be defined in `mesh.inp` for the volume computation to work — it represents the lumen-facing nodes within the fibrotic seed region.
- All scripts assume the working directory contains the relevant `.odb`, `.inp`, and `.csv` files.
