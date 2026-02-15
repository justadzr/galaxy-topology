import numpy as np
import matplotlib.pyplot as plt
import os

# Fix font issues
# Using internal mathtext parser (not external TeX) to avoid finding executables
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'serif' # Use a serif font for math-like appearance
plt.rcParams['mathtext.fontset'] = 'cm' # Computer Modern for math
plt.rcParams['axes.unicode_minus'] = False # Fix minus sign issues


def plot_betti_curves(npz_file, output_file):
    if not os.path.exists(npz_file):
        print(f"Error: {npz_file} not found.")
        return

    data = np.load(npz_file)
    print(f"Keys in npz: {list(data.keys())}")
    
    # bc_scaled shape: (N_selections, 3, resolution)
    bc_scaled = data['bc_scaled']
    
    # Try different ways to get alpha
    if 'alpha_scaled' in data:
        alpha_scaled = data['alpha_scaled']
    else:
        # Assuming 0 to 2 range and shape of bc_scaled
        resolution = bc_scaled.shape[-1]
        alpha_scaled = np.linspace(0, 2, resolution)
        
    print(f"DEBUG: alpha_scaled shape: {alpha_scaled.shape}")
    print(f"DEBUG: bc_scaled shape: {bc_scaled.shape}")
    
    # Handle extra dimension (list of simulations)
    if bc_scaled.ndim == 4 and bc_scaled.shape[0] == 1:
         bc_scaled = bc_scaled[0]
         print(f"DEBUG: Squeezed bc_scaled shape: {bc_scaled.shape}")
    
    # N_selections for TNG usually 4 (different number densities)
    n_selections = bc_scaled.shape[0]
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    dims = [0, 1, 2]
    colors = ['r', 'g', 'b', 'k'] # One for each selection
    linestyles = ['-', '--', '-.', ':']
    # Corrected units: density is usually in (Mpc/h)^-3 = h^3 Mpc^-3? 
    # Or number density n is in (h^-1 Mpc)^-3 which is h^3 Mpc^-3.
    # The previous label was "1e-2 (h/Mpc)^3" which is h^3 Mpc^-3.
    # The user suggests "h^-1 Mpc" for length, so density is (h^-1 Mpc)^-3.
    # This is equivalent to h^3 Mpc^-3.
    # Let's use the explicit sample names to match the paper
    labels = ["Halos", "All Galaxies", "Star-Forming", "Quiescent"]

    # Plot each dimension in a subplot
    for d_idx, d in enumerate(dims):
        ax = axes[d_idx]
        for i in range(n_selections):
            if i >= len(colors): break # Safety
            
            # bc_scaled[i, d, :] is the curve
            label = labels[i] if i < len(labels) else f"Sel {i}"
            
            # Important: check if axis order is correct in loaded npz
            # Shape is (N_sel, 3, Res). Check dim d.
            
            curve = bc_scaled[i, d, :]
            
            # Ensure 1D
            x = alpha_scaled.flatten()
            y = curve.flatten()
            
            if x.shape[0] != y.shape[0]:
                print(f"Error: Shape mismatch for Sel {i} Dim {d}: x={x.shape}, y={y.shape}")
                continue
                
            ax.plot(x, y, color=colors[i], linestyle=linestyles[i], label=label)
        
        ax.set_title(f"$\\beta_{d}$")
        # ax.set_xlabel(r"$\nu = \alpha \bar{n}^{1/3}$") # Original confusing label
        ax.set_xlabel(r"$\nu = \alpha n^{1/3}$") # Removed bar to avoid confusion with tilde
        ax.set_xlim(0, 2)
        ax.grid(True, alpha=0.3)
        if d_idx == 0:
            ax.set_ylabel(r"$\beta_k / N$") # Normalized per object
            ax.legend(loc="upper right")

    plt.suptitle(f"TNG-100 z=0 Topology (File: es_all_tng_z0_final_v13.npz)")
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    npz_path = "/homes/yzang27/tda/galaxy-topology/topology_summaries/TNG-100/es_all_ssfr_match_snap91.npz"
    output_path = "tng_topology_ssfr_match_snap91.png"
    plot_betti_curves(npz_path, output_path)
    # data = np.load("/homes/yzang27/tda/galaxy-topology/topology_summaries/TNG-100/es_all_ssfr_match.npz")
    # print(data['bc_scaled'].shape)