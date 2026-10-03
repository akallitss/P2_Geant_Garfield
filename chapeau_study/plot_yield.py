#!/usr/bin/env python3
"""Photon yield vs photon energy leaving an Al chapeau wall (119 MeV e-, 10 deg).
Nominal wall 2 mm (see wall_estimate.py), 1-5 mm band.  Yield = dN/dk per
incident electron per keV."""
import math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import chapeau_photons as cp

E, INC = 119.0, 10
k = np.geomspace(3, 150, 400)
s = lambda t_mm: t_mm/10/math.cos(math.radians(INC))
per_keV = lambda dN_dlnk: dN_dlnk / k
C = {"1": "#1b9e77", "2": "#d95f02", "5": "#7570b3"}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

# ---- plot 1: total exiting yield, thickness band ----------------------------
fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
for t in (1, 2, 5):
    tot = cp.primary_spectrum(k, E, s(t)) + cp.delta_spectrum(k, E, s(t))
    ax.loglog(k, per_keV(tot), c=C[str(t)], lw=2.5 if t == 2 else 1.5,
              label=f"{t} mm Al" + ("  (nominal)" if t == 2 else ""))
ax.set(xlabel="photon energy [keV]", ylabel="photons / electron / keV",
       title="Photons exiting the Al wall, 119 MeV e⁻ at 10°", xlim=(3, 150), ylim=(1e-6, 1e-1))
ax.grid(alpha=.3, which="both"); ax.legend(frameon=False)
fig.savefig("yield_vs_energy.png", dpi=160)

# ---- plot 2: by process, produced vs exiting (2 mm) ---------------------------
t = 2; ss = s(t)
prim_exit = cp.primary_spectrum(k, E, ss)
delt_exit = cp.delta_spectrum(k, E, ss)
# "produced": same formulas with the wall made transparent
orig = cp.mu_al; cp.mu_al = lambda kk: np.full_like(np.asarray(kk, float), 1e-9)
prim_prod = cp.primary_spectrum(k, E, ss)
delt_prod = cp.delta_spectrum(k, E, ss)
cp.mu_al = orig
fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
ax.loglog(k, per_keV(prim_prod), c="#1b9e77", ls=":", label="primary brems, produced")
ax.loglog(k, per_keV(prim_exit), c="#1b9e77", lw=2, label="primary brems, exits wall")
ax.loglog(k, per_keV(delt_prod), c="#d95f02", ls=":", label="δ-ray brems, produced")
ax.loglog(k, per_keV(delt_exit), c="#d95f02", lw=2, label="δ-ray brems, exits wall")
ax.set(xlabel="photon energy [keV]", ylabel="photons / electron / keV",
       title="By process, 2 mm Al, 119 MeV e⁻ at 10°", xlim=(3, 150), ylim=(1e-6, 1e-1))
ax.grid(alpha=.3, which="both"); ax.legend(frameon=False, fontsize=8)
fig.savefig("yield_by_process.png", dpi=160)
print("done")
