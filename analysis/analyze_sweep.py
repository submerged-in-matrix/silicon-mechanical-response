import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

def read_ss(filename):
    """Read stress-strain data file."""
    data = []
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    data.append([float(parts[1]), float(parts[2])])
                except ValueError:
                    continue
    return np.array(data)

def get_E(strain, stress, lo=0.03, hi=0.08):
    """Extract Young's modulus from a strain region."""
    w = 30
    if len(stress) < w + 10:
        w = max(3, len(stress) // 5)
    stress_s = np.convolve(stress, np.ones(w)/w, mode='valid')
    strain_s = strain[(w-1)//2 : (w-1)//2 + len(stress_s)]
    mask = (strain_s > lo) & (strain_s < hi)
    if mask.sum() < 5:
        return 0, 0, strain_s, stress_s
    c = np.polyfit(strain_s[mask], stress_s[mask], 1)
    pred = np.polyval(c, strain_s[mask])
    ss_res = np.sum((stress_s[mask] - pred) ** 2)
    ss_tot = np.sum((stress_s[mask] - stress_s[mask].mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    return c[0], r2, strain_s, stress_s

# --- Parameters ---
radii = [3, 5, 7, 10]
rates = [0.0001, 0.0005, 0.001, 0.005]
lattice_a = 5.431  # Angstrom

# Atom counts from the runs
atom_counts = {3: 2700, 5: 7596, 7: 14748, 10: 30252}

# --- Scan best fit region for each run ---
fit_regions = [
    (0.01, 0.03), (0.02, 0.05), (0.03, 0.06),
    (0.03, 0.08), (0.04, 0.08), (0.05, 0.10),
]

# --- Extract E for all combinations ---
results = {}  # (rad, rate) -> E

print(f"{'RAD':>4} {'Diameter(A)':>12} {'Atoms':>7} {'Rate':>10} {'E(GPa)':>8} {'R2':>6} {'Region':>8}")
print("-" * 65)

for rad in radii:
    for rate in rates:
        fname = f'../outputs/ss_R{rad}_E{rate}.dat'
        if not os.path.exists(fname):
            continue
        data = read_ss(fname)
        if len(data) < 50:
            continue

        strain = data[:, 0]
        stress = data[:, 1]

        # Find best fit region
        best_E = 0
        best_r2 = 0
        best_reg = ""
        for lo, hi in fit_regions:
            E, r2, _, _ = get_E(strain, stress, lo, hi)
            if r2 > 0.98 and E > best_E:
                best_E = E
                best_r2 = r2
                best_reg = f"{lo*100:.0f}-{hi*100:.0f}%"

        if best_E == 0:
            for lo, hi in fit_regions:
                E, r2, _, _ = get_E(strain, stress, lo, hi)
                if r2 > best_r2:
                    best_E = E
                    best_r2 = r2
                    best_reg = f"{lo*100:.0f}-{hi*100:.0f}%"

        diameter = 2 * rad * lattice_a
        results[(rad, rate)] = best_E
        print(f"{rad:4d} {diameter:12.1f} {atom_counts[rad]:7d} {rate:10.4f} {best_E:8.1f} {best_r2:6.3f} {best_reg:>8}")

# --- Summary tables ---
print(f"\n{'='*60}")
print("YOUNG'S MODULUS (GPa) - RADIUS vs STRAIN RATE")
print(f"{'='*60}")
print(f"{'':>12}", end="")
for rate in rates:
    print(f"  {rate:>10}", end="")
print(f"  {'Diameter':>10}")
print("-" * 70)
for rad in radii:
    diameter = 2 * rad * lattice_a
    print(f"RAD={rad:>2}   ", end="")
    for rate in rates:
        E = results.get((rad, rate), 0)
        print(f"  {E:10.1f}", end="")
    print(f"  {diameter:8.1f} A")
print(f"\nSi [100] bulk experimental: ~130 GPa")
print(f"{'='*60}")

# --- Figure 1: E vs diameter for each strain rate ---
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

ax = axes[0]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
for i, rate in enumerate(rates):
    diameters = [2 * r * lattice_a for r in radii]
    E_vals = [results.get((r, rate), 0) for r in radii]
    ax.plot(diameters, E_vals, 'o-', color=colors[i], markersize=10,
            linewidth=2, label=f'rate = {rate}')

ax.axhline(y=130, color='black', linestyle=':', linewidth=1.5,
           label='Bulk Si [100] = 130 GPa')
ax.set_xlabel('Wire Diameter ($\\AA$)', fontsize=14)
ax.set_ylabel("Young's Modulus (GPa)", fontsize=14)
ax.set_title('(a) Size Effect', fontsize=14)
ax.legend(fontsize=9)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

# --- Figure 2: E vs strain rate for each radius ---
ax = axes[1]
for i, rad in enumerate(radii):
    diameter = 2 * rad * lattice_a
    rate_vals = rates
    E_vals = [results.get((rad, r), 0) for r in rates]
    ax.semilogx(rate_vals, E_vals, 's-', color=colors[i], markersize=10,
                linewidth=2, label=f'd = {diameter:.0f} $\\AA$')

ax.axhline(y=130, color='black', linestyle=':', linewidth=1.5,
           label='Bulk Si [100]')
ax.set_xlabel('Strain Rate (1/ps)', fontsize=14)
ax.set_ylabel("Young's Modulus (GPa)", fontsize=14)
ax.set_title('(b) Strain Rate Effect', fontsize=14)
ax.legend(fontsize=9)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

# --- Figure 3: Stress-strain curves for largest wire at all rates ---
ax = axes[2]
for i, rate in enumerate(rates):
    fname = f'../outputs/ss_R10_E{rate}.dat'
    if not os.path.exists(fname):
        continue
    data = read_ss(fname)
    strain = data[:, 0]
    stress = data[:, 1]
    w = 30
    stress_s = np.convolve(stress, np.ones(w)/w, mode='valid')
    strain_s = strain[(w-1)//2 : (w-1)//2 + len(stress_s)]
    ax.plot(strain_s * 100, stress_s, '-', color=colors[i], linewidth=1.5,
            label=f'rate = {rate}')

ax.set_xlabel('Engineering Strain (%)', fontsize=14)
ax.set_ylabel('Stress (GPa)', fontsize=14)
ax.set_title('(c) Stress-Strain: R=10 (d=108$\\AA$)', fontsize=14)
ax.legend(fontsize=10)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
os.makedirs('../figures', exist_ok=True)
plt.savefig('../figures/parametric_sweep.png', dpi=300, bbox_inches='tight')
print(f"\nFigure saved: ../figures/parametric_sweep.png")