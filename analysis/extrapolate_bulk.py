import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# --- Data from our sweep (averaged over strain rates) ---
lattice_a = 5.431  # Angstrom
radii = [3, 5, 7, 10]  # in lattice units

# E values from sweep (all 4 rates per radius)
E_all = {
    3: [35.9, 33.8, 32.7, 34.6],
    5: [37.0, 38.5, 37.7, 40.7],
    7: [39.3, 40.6, 41.0, 42.4],
    10: [42.8, 43.2, 43.1, 44.5],
}

# Average E and std error for each radius
diameters = []
E_mean = []
E_err = []

for rad in radii:
    d = 2 * rad * lattice_a
    diameters.append(d)
    E_mean.append(np.mean(E_all[rad]))
    E_err.append(np.std(E_all[rad]) / np.sqrt(len(E_all[rad])))

diameters = np.array(diameters)
E_mean = np.array(E_mean)
E_err = np.array(E_err)

# --- Surface fraction estimate ---
# For a cylinder of radius R, surface shell thickness ~ one bond length
# In diamond cubic, nearest neighbor distance = a*sqrt(3)/4 = 2.35 A
shell = lattice_a  # conservative: one lattice constant thick
radii_angstrom = diameters / 2.0
surf_frac = 1.0 - ((radii_angstrom - shell) / radii_angstrom) ** 2
# Clamp negative values for very small wires
surf_frac = np.clip(surf_frac, 0, 1)

print("=" * 65)
print("SURFACE FRACTION AND BULK EXTRAPOLATION")
print("=" * 65)
print(f"{'Diameter(A)':>12} {'Atoms':>8} {'E_avg(GPa)':>10} {'E_err':>8} {'Surf_frac':>10}")
print("-" * 65)
atom_counts = [2700, 7596, 14748, 30252]
for i in range(len(diameters)):
    print(f"{diameters[i]:12.1f} {atom_counts[i]:8d} {E_mean[i]:10.1f} {E_err[i]:8.2f} {surf_frac[i]:10.1%}")

# --- Extrapolation 1: E vs 1/diameter ---
inv_d = 1.0 / diameters
coeffs_d = np.polyfit(inv_d, E_mean, 1)
E_bulk_from_d = coeffs_d[1]  # intercept = E at 1/d -> 0

print(f"\nLinear fit: E = {coeffs_d[1]:.1f} + {coeffs_d[0]:.1f} / d")
print(f"Extrapolated bulk E (1/d -> 0): {E_bulk_from_d:.1f} GPa")

# --- Extrapolation 2: E vs surface fraction ---
coeffs_sf = np.polyfit(surf_frac, E_mean, 1)
E_bulk_from_sf = coeffs_sf[1]  # intercept = E at surf_frac -> 0

print(f"\nLinear fit: E = {coeffs_sf[1]:.1f} + {coeffs_sf[0]:.1f} * surf_frac")
print(f"Extrapolated bulk E (surf -> 0): {E_bulk_from_sf:.1f} GPa")
print(f"Experimental Si [100]:           130 GPa")
print("=" * 65)

# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: E vs 1/diameter
ax = axes[0]
ax.errorbar(inv_d * 1000, E_mean, yerr=E_err, fmt='ro', markersize=10,
            capsize=5, linewidth=2, label='MD data')

# Fit line extended to 1/d = 0
inv_d_fit = np.linspace(0, max(inv_d) * 1.1, 100)
E_fit = np.polyval(coeffs_d, inv_d_fit)
ax.plot(inv_d_fit * 1000, E_fit, 'b--', linewidth=2,
        label=f'Linear fit (E$_{{bulk}}$ = {E_bulk_from_d:.0f} GPa)')

ax.axhline(y=130, color='green', linestyle=':', linewidth=2,
           label='Expt bulk Si [100] = 130 GPa')
ax.plot(0, E_bulk_from_d, 'b*', markersize=20, label='Extrapolated bulk')

ax.set_xlabel('1 / Diameter (1000 / $\\AA$)', fontsize=14)
ax.set_ylabel("Young's Modulus (GPa)", fontsize=14)
ax.set_title('(a) Bulk Extrapolation via 1/d', fontsize=14)
ax.legend(fontsize=10)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

# Right: E vs surface fraction
ax = axes[1]
ax.errorbar(surf_frac * 100, E_mean, yerr=E_err, fmt='rs', markersize=10,
            capsize=5, linewidth=2, label='MD data')

sf_fit = np.linspace(0, max(surf_frac) * 1.1, 100)
E_sf_fit = np.polyval(coeffs_sf, sf_fit)
ax.plot(sf_fit * 100, E_sf_fit, 'b--', linewidth=2,
        label=f'Linear fit (E$_{{bulk}}$ = {E_bulk_from_sf:.0f} GPa)')

ax.axhline(y=130, color='green', linestyle=':', linewidth=2,
           label='Expt bulk Si [100] = 130 GPa')
ax.plot(0, E_bulk_from_sf, 'b*', markersize=20, label='Extrapolated bulk')

ax.set_xlabel('Estimated Surface Fraction (%)', fontsize=14)
ax.set_ylabel("Young's Modulus (GPa)", fontsize=14)
ax.set_title('(b) E vs Surface Atom Fraction', fontsize=14)
ax.legend(fontsize=10)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
os.makedirs('../figures', exist_ok=True)
plt.savefig('../figures/bulk_extrapolation.png', dpi=300, bbox_inches='tight')
print(f"\nFigure saved: ../figures/bulk_extrapolation.png")
