import numpy as np
import h5py
import gudhi as gd
import gudhi.representations as gdr
# import homcloud.interface as hc
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.colors import LogNorm
from scipy import stats
from scipy.interpolate import interp1d
from Corrfunc.theory.xi import xi
from alpha_complex_periodic import calc_persistence
import illustris_python as tng

mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})

snap = 91

with np.load(f"/homes/yzang27/tda/galaxy-topology/topology_summaries/TNG-100/es_all_snap{snap}.npz") as data:
    alpha_tng300 = data["alpha"]
    alpha_scaled = data["alpha_scaled"]
    bc_tng300 = data["bc"]
    bc_scaled_tng300 = data["bc_scaled"]
    xi_ravg_tng300 = data["xi_ravg"]
    xi_tng300 = data["xi"]

ylims1 = [5e-4, 5e-5, 7e-6]
ylims2 = [2e-3, 5e-4, 6e-5]

fig, axs = plt.subplots(3, 2, figsize=(3.38,4), sharex="col", constrained_layout=True, dpi=200)
for d in range(3):

    ax = axs[d,0]
    ax.plot(alpha_tng300, bc_tng300[0,1,d,:], 'b', alpha=1, label="TNG100")
    ax.loglog()
    ax.set_xlim(2e-2, 15)
    ax.set_ylim(bottom=ylims1[d])
    if not d: ax.legend(fontsize=6, frameon=False)
    ax.set_ylabel(f"$\\beta_{d}$ [$h^3$ Mpc$^{{-3}}$]")
    
    ax = axs[d,1]
    ax.plot(alpha_scaled, bc_scaled_tng300[0,1,d,:], 'b', alpha=1)
    ax.loglog()
    ax.set_xlim(6e-3, 5)
    ax.set_ylim(bottom=ylims2[d])
    ax.set_ylabel(f"$\\beta_{d}$ [$\ell^{{-3}}$]")
    
axs[2,0].set_xlabel("$\\alpha$ [$h^{-1}$ Mpc]")
axs[2,1].set_xlabel("$\\alpha$ [$\ell$]")

plt.savefig(f"/homes/yzang27/tda/galaxy-topology/pictures/TNG-100/galaxies_betti_full_vert_{snap}.png")
# plt.show()

xlims = [1, 2.5, 4.5]
ylims = [2e-3, 5e-4, 5e-5]

def mean_conf_int(data, p=0.68):
    mean = np.mean(data, axis=0)
    lq = np.quantile(data, 0.5 - p/2, axis=0)
    hq = np.quantile(data, 0.5 + p/2, axis=0)
    return mean, lq, hq

fig, axs = plt.subplots(3, 1, figsize=(2,4), sharex=True, constrained_layout=True, dpi=200)
for d in range(3):
    ax = axs[d]
    # mean_gal, lq_gal, hq_gal = mean_conf_int(bc_scaled_cv_sam[:,1,d,:], p=0.95)
    # mean_sf, lq_sf, hq_sf = mean_conf_int(bc_scaled_cv_sam[:,2,d,:], p=0.95)
    # mean_qsnt, lq_qsnt, hq_qsnt = mean_conf_int(bc_scaled_cv_sam[:,3,d,:], p=0.95)
    
    mean_gal_tng300, lq_gal_tng300, hq_gal_tng300 = mean_conf_int(bc_scaled_tng300[:,1,d,:], p=0.95)
    mean_sf_tng300, lq_sf_tng300, hq_sf_tng300 = mean_conf_int(bc_scaled_tng300[:,2,d,:], p=0.95)
    mean_qsnt_tng300, lq_qsnt_tng300, hq_qsnt_tng300 = mean_conf_int(bc_scaled_tng300[:,3,d,:], p=0.95)

    # ax.plot(alpha_scaled, mean_gal, color='blue', label="all galaxies" if not d else "CAMELS-SAM")
    # ax.fill_between(alpha_scaled, lq_gal, hq_gal, color='blue', alpha=0.2)
    # ax.plot(alpha_scaled, mean_sf, color='orange', label="star-forming" if not d else "")
    # ax.fill_between(alpha_scaled, lq_sf, hq_sf, alpha=0.3, color='orange')
    # ax.plot(alpha_scaled, mean_qsnt, color='green', label="quiescent" if not d else "")
    # ax.fill_between(alpha_scaled, lq_qsnt, hq_qsnt, color='green', alpha=0.3)
    
    ax.plot(alpha_scaled, mean_gal_tng300, color='blue', label="all galaxies" if not d else "")
    ax.plot(alpha_scaled, mean_sf_tng300, color='orange', label="star-forming" if not d else "")
    ax.plot(alpha_scaled, mean_qsnt_tng300, color='green', label="quiescent" if not d else "")

    ax.loglog()
    ax.set_xlim(1e-2, 5)
    #ticks = [0.01, 0.1, 1, 2, 4]
    #plt.xticks(ticks, ticks)
    ax.set_ylim(bottom=ylims[d])
    if not d: ax.legend(fontsize=6, frameon=False)
    #if d == 1: plt.legend(fontsize=6, frameon=False)

    ax.set_ylabel(f"$\\beta_{d}$ [$\ell^{{-3}}$]")
axs[2].set_xlabel("$\\alpha$ [$\\ell$]")
    
plt.savefig(f"/homes/yzang27/tda/galaxy-topology/pictures/TNG-100/galaxies_betti_vert_{snap}.png")
# plt.show()