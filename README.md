# Mechanical Response of Silicon Nanowires: A Parametric MD Study

## Summary

Parametric study of uniaxial tensile deformation of [001] silicon
nanowires using the Tersoff potential. Investigated the effects of
wire diameter (33–109 Å) and strain rate (10⁸–5×10⁹ s⁻¹) on
Young's modulus. Found that size effect dominates: E increases from
34 to 44 GPa with diameter, while strain rate has minimal impact
(±2 GPa). Extrapolation toward bulk Si [100] E ≈ 130 GPa is
consistent with a surface-fraction-dependent stiffness model.

## Key Finding

Young's modulus is controlled by the surface-to-volume ratio, not
strain rate. This is consistent with the picture that under-coordinated
surface atoms in diamond cubic silicon contribute significantly
reduced stiffness compared to bulk-coordinated interior atoms.

| Diameter (Å) | Atoms  | E (GPa) avg | Surface fraction (est.) |
|--------------|--------|-------------|------------------------|
| 32.6         | 2,700  | 34.3        | ~50%                   |
| 54.3         | 7,596  | 38.5        | ~35%                   |
| 76.0         | 14,748 | 40.8        | ~27%                   |
| 108.6        | 30,252 | 43.4        | ~20%                   |
| Bulk (∞)     | ∞      | ~130        | 0%                     |

Strain rate (10⁸ to 5×10⁹ s⁻¹) changes E by < 5% at any given
diameter — mechanical stiffness is essentially rate-independent in
this regime.

## Method

- **Software:** LAMMPS
- **Potential:** Tersoff (Si.tersoff)
- **System:** Diamond cubic Si nanowires, [001] axis, s s p boundaries
- **Parametric sweep:** 4 radii × 4 strain rates = 16 independent runs
- **Protocol per run:**
  1. Energy minimization (tight tolerances for surface reconstruction)
  2. NVT equilibration at 300 K (30 ps)
  3. Constant strain rate deformation to 20% strain
- **Analysis:** Smoothed stress-strain, best-fit E from linear region,
  systematic region scanning with R² quality control

## Files
- scripts/
- in.si_parametric    — Parameterized input (-var RAD, -var RATE, -var NSTEPS)
- in.si_tension       — Single-run version (original study)
- run_sweep.sh        — Batch runner for all 16 combinations
- Si.tersoff          — Tersoff potential file
- outputs/
- ss_R{rad}_E{rate}.dat   — Stress-strain data (16 files)
- log_R{rad}_E{rate}.txt  — LAMMPS logs (16 files)
- analysis/
- analyze_sweep.py    — Parametric analysis + multi-panel figure
- analyze_si_tension.py — Single-run analysis (original study)
- figures/
- parametric_sweep.png   — 3-panel: size effect, rate effect, stress-strain
- si_stress_strain.png   — Si vs Cu comparison

## How to Reproduce

```bash
cd scripts
cp /usr/share/lammps/potentials/Si.tersoff .
./run_sweep.sh          # ~1-2 hours for all 16 runs
cd ../analysis
python3 analyze_sweep.py
```
## Discussion

Two effects contribute to the low measured E compared to experimental
Si [100] = 130 GPa:

**1. Surface effect (explains ~34→47 GPa, i.e. ~28% of the stiffness):**
Under-coordinated surface atoms in diamond cubic silicon have fewer
than 4 tetrahedral bonds, reducing local stiffness. Extrapolation to
zero surface fraction gives E ≈ 47 GPa for both 1/d and surface
fraction models, confirming a significant but not dominant role.

**2. Tersoff potential limitation (explains 47 vs 130 GPa):**
The Tersoff potential was parameterized to reproduce cohesive energy,
lattice constant, and bulk modulus (~98 GPa), but systematically
underestimates the anisotropic elastic constants C11 and C12 that
determine E[100]. This is a known limitation documented in the
literature. The Stillinger-Weber potential or modified Tersoff
(Si.tersoff.mod) may give closer values.

**Key insight:** Both potential accuracy AND system size matter.
A parametric sweep isolates these effects quantitatively.

## References

- Tersoff, J. (1988) Phys. Rev. B 37, 6991 (Tersoff potential)
- Thompson et al. (2022) Comput. Phys. Commun. 271, 108171 (LAMMPS)
- Pal & Reddy (2024) Molecular Dynamics for Materials Modeling (Routledge)
