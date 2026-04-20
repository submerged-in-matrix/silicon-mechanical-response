import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# Read data
data = []
with open('../outputs/stress_strain.dat') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                data.append([float(parts[1]), float(parts[2])])
            except ValueError:
                continue

data = np.array(data)
strain = data[:, 0]
stress = data[:, 1]

# Smooth
w = 30
stress_s = np.convolve(stress, np.ones(w)/w, mode='valid')
strain_s = strain[(w-1)//2 : (w-1)//2 + len(stress_s)]

# Scan strain regions for E
print(f"{'Region':>16} {'E (GPa)':>10} {'R-squared':>10}")
print("-" * 40)

regions = [
    (0.005, 0.02, "0.5-2%"),
    (0.01,  0.03, "1-3%"),
    (0.02,  0.04, "2-4%"),
    (0.03,  0.05, "3-5%"),
    (0.04,  0.06, "4-6%"),
    (0.05,  0.08, "5-8%"),
    (0.02,  0.06, "2-6%"),
    (0.03,  0.08, "3-8%"),
]

best_E = 0
best_region = ""
best_r2 = 0
best_coeffs = None

for lo, hi, name in regions:
    mask = (strain_s > lo) & (strain_s < hi)
    if mask.sum() < 5:
        continue
    c = np.polyfit(strain_s[mask], stress_s[mask], 1)
    E = c[0]
    pred = np.polyval(c, strain_s[mask])
    ss_res = np.sum((stress_s[mask] - pred) ** 2)
    ss_tot = np.sum((stress_s[mask] - stress_s[mask].mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    print(f"{name:>16} {E:10.1f} {r2:10.4f}")
    if r2 > 0.99 and E > best_E:
        best_E = E
        best_region = name
        best_r2 = r2
        best_coeffs = c

if best_coeffs is None:
    for lo, hi, name in regions:
        mask = (strain_s > lo) & (strain_s < hi)
        if mask.sum() < 5:
            continue
        c = np.polyfit(strain_s[mask], stress_s[mask], 1)
        pred = np.polyval(c, strain_s[mask])
        ss_res = np.sum((stress_s[mask] - pred) ** 2)
        ss_tot = np.sum((stress_s[mask] - stress_s[mask].mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        if r2 > best_r2:
            best_r2 = r2
            best_E = c[0]
            best_region = name
            best_coeffs = c

# Peak stress
peak_idx = np.argmax(stress_s)
peak_strain = strain_s[peak_idx]
peak_stress = stress_s[peak_idx]

print(f"\nBest fit region: {best_region}")
print(f"\n{'='*55}")
print(f"SILICON NANOWIRE TENSILE TEST (Tersoff)")
print(f"{'='*55}")
print(f"System: 4050 atoms, diamond cubic, [001] axis")
print(f"Strain rate: 5 x 10^8 s^-1, Temperature: 300 K")
print(f"{'-'*55}")
print(f"Young's modulus:     {best_E:.1f} GPa ({best_region})")
print(f"Si [100] expt:       ~130 GPa")
print(f"Peak stress:         {peak_stress:.1f} GPa at {peak_strain*100:.1f}%")
print(f"Si bulk fracture:    ~7 GPa (theoretical)")
print(f"{'='*55}")

# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Left: full curve
ax = axes[0]
ax.plot(strain * 100, stress, 'b-', linewidth=0.6, alpha=0.3, label='Raw')
ax.plot(strain_s * 100, stress_s, 'b-', linewidth=2, label='Smoothed')

for lo, hi, name in regions:
    if name == best_region:
        fit_x = np.linspace(0, hi + 0.02, 100)
        fit_y = np.polyval(best_coeffs, fit_x)
        ax.plot(fit_x * 100, fit_y, 'r--', linewidth=2.5,
                label=f'E = {best_E:.0f} GPa ({best_region})')
        break

ax.plot(peak_strain * 100, peak_stress, 'rv', markersize=14,
        label=f'Peak: {peak_stress:.1f} GPa')

ax.set_xlabel('Engineering Strain (%)', fontsize=14)
ax.set_ylabel('Stress (GPa)', fontsize=14)
ax.set_title('(a) Si Nanowire [001] Tension', fontsize=14)
ax.legend(fontsize=10)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

# Right: comparison with Cu nanowire
ax = axes[1]
ax.plot(strain_s * 100, stress_s, 'b-', linewidth=2, label='Si (Tersoff, this work)')

# Load Cu data if available
cu_file = '../../MD_LAMMPS_Ex/MD_LAMMPS_Ex/project-02-nanowire-tension/outputs/stress_strain.dat'
if os.path.exists(cu_file):
    cu_data = []
    with open(cu_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    cu_data.append([float(parts[1]), float(parts[2])])
                except ValueError:
                    continue
    cu_data = np.array(cu_data)
    cu_strain = cu_data[:, 0]
    cu_stress = cu_data[:, 1]
    cu_w = 30
    cu_stress_s = np.convolve(cu_stress, np.ones(cu_w)/cu_w, mode='valid')
    cu_strain_s = cu_strain[(cu_w-1)//2 : (cu_w-1)//2 + len(cu_stress_s)]
    ax.plot(cu_strain_s * 100, cu_stress_s, 'r-', linewidth=2,
            label='Cu (EAM, Project 2)')

ax.set_xlabel('Engineering Strain (%)', fontsize=14)
ax.set_ylabel('Stress (GPa)', fontsize=14)
ax.set_title('(b) Si vs Cu: Brittle vs Ductile', fontsize=14)
ax.legend(fontsize=11)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
os.makedirs('../figures', exist_ok=True)
plt.savefig('../figures/si_stress_strain.png', dpi=300, bbox_inches='tight')
print(f"\nFigure saved")
