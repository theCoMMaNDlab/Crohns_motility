## UMAT parameter correspondence

The UMAT always expects exactly 3 constants, in this fixed order, on the
`*USER MATERIAL` line of the `.inp` file:

```
PROPS(1) = KF_CIRC   circumferential fiber stiffness
PROPS(2) = KF_LONG   longitudinal fiber stiffness
PROPS(3) = MU0       isotropic (neo-Hookean) shear modulus
```

## Commands

### Circumferential

```
abaqus job=single_elem_circum user=umat_kcalib.f interactive
abaqus python extract_and_plot_stress_strain.py --extract single_elem_circum circum
python extract_and_plot_stress_strain.py --score circum_job_response.txt circum_stretch_stress.txt circum
```

### Longitudinal

```
abaqus job=single_elem_long user=umat_kcalib.f interactive
abaqus python extract_and_plot_stress_strain.py --extract single_elem_long long
python extract_and_plot_stress_strain.py --score long_job_response.txt long_stretch_stress.txt long
```
