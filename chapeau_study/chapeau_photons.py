#!/usr/bin/env python3
"""
Analytic estimate: can ~120-145 MeV back-scattered electrons crossing an
aluminium "chapeau" wall produce the low-energy photons seen in simulation?

No Geant4 here -- this is the pencil-and-paper expectation a real simulation
should be compared against.  Photon sources inside the Al slab:

  (1) primary bremsstrahlung of the beam-like electron   (Bethe-Heitler, 1/k)
  (2) bremsstrahlung of knock-on (delta) electrons       (soft, many of them)

Not modelled (all expected to be minor in a thin Al wall): Compton/photoeff.
of the hard brems photons, Al K-fluorescence (1.49 keV, dies in microns),
transition radiation, annihilation, backscatter from things behind the wall.

Inputs that are GUESSES: wall thickness (scanned), incidence angle (scanned).
Data: NIST XCOM mass attenuation (fetched 2026-10-03, physics.nist.gov
XrayMassCoef z13 / z18); X0(Al) = 24.01 g/cm2 (PDG).

Usage:  python chapeau_photons.py        (writes table + figures here)
"""
import math
from pathlib import Path
import numpy as np

OUT = Path(__file__).parent

# ---- constants / materials -------------------------------------------------
ME = 0.511                      # MeV
RHO_AL, X0_AL_G = 2.699, 24.01  # g/cm3, g/cm2
X0_AL = X0_AL_G / RHO_AL        # cm  (8.9 cm)
ZA_AL = 13 / 26.9815
K_DELTA = 0.1535 * ZA_AL * RHO_AL   # MeV/cm : dN/dT = K/(beta^2 T^2) per cm
RHO_AR = 1.66e-3                # g/cm3, 1 atm 20C (pure Ar; mixture unknown)

# NIST: E [keV], mu/rho [cm2/g]
AL = np.array([
 (1,1185),(1.5,402.2),(2,2263),(3,788.0),(4,360.5),(5,193.4),(6,115.3),
 (8,50.33),(10,26.23),(15,7.955),(20,3.441),(30,1.128),(40,0.5685),
 (50,0.3681),(60,0.2778),(80,0.2018),(100,0.1704),(150,0.1378),(200,0.1223)])
# NIST Ar: E [keV], mu_en/rho [cm2/g]  (energy-transfer weighting)
AR_EN = np.array([
 (2,509.3),(3,168.2),(4,697.9),(5,395.3),(6,244.9),(8,112.5),(10,60.38),
 (15,18.86),(20,8.074),(30,2.382),(40,0.9907),(50,0.5020),(60,0.2904),
 (80,0.1280),(100,0.07344),(150,0.03703),(200,0.02998)])

def loglog(tab, k):
    return np.exp(np.interp(np.log(k), np.log(tab[:, 0]), np.log(tab[:, 1])))

def mu_al(k_kev):   # 1/cm
    return loglog(AL, k_kev) * RHO_AL

# ---- electron range (Katz-Penfold, MeV) -----------------------------------
def csda_range_cm(T):           # kinetic T in MeV, Al
    T = np.clip(T, 0.01, 3.0)
    return 0.412 * T**(1.265 - 0.0954*np.log(T)) / RHO_AL

# ---- bremsstrahlung shape --------------------------------------------------
def f_bh(y):
    """Bethe-Heitler (complete screening, no LPM) shape: 4/3(1-y)+y^2."""
    return np.where(y < 1, 4/3*(1-y) + y*y, 0.0)

# ---- (1) primary brems -----------------------------------------------------
def primary_spectrum(k_kev, E_MeV, s_cm):
    """dN/dlnk per electron escaping a slab of path s.  Emission uniform in
    depth, photons travel ~parallel to the electron (cone ~ m/E << MS angle)."""
    k = np.asarray(k_kev) / 1000.0
    mu = mu_al(k_kev)
    esc = (1 - np.exp(-mu*s_cm)) / mu          # cm, = integral of exp(-mu(s-x))dx
    return f_bh(k/(E_MeV+ME)) * esc / X0_AL

# ---- (2) delta-ray brems -----------------------------------------------------
def delta_spectrum(k_kev, E_MeV, s_cm, Tmin=0.010, nT=60):
    """dN/dlnk per primary from brems of knock-on electrons (T>Tmin).
    Each delta radiates over min(range, s/2) with energy loss followed;
    photons attenuated over s/2 on average.  Bethe-Heitler shape used down to
    ~10 keV where it is only good to a factor ~2 (hence 'estimate')."""
    gam = (E_MeV+ME)/ME; beta2 = 1 - 1/gam**2
    Tmax = E_MeV / 2
    Ts = np.geomspace(Tmin, Tmax, nT)
    dT = np.gradient(Ts)
    k = np.asarray(k_kev)/1000.0
    out = np.zeros_like(k)
    att = np.exp(-mu_al(k_kev) * s_cm/2)
    for T0, w in zip(Ts, dT):
        n_delta = K_DELTA / (beta2*T0**2) * s_cm * w     # deltas created
        cap = min(float(csda_range_cm(T0)), s_cm/2)
        # follow slowing down in 20 steps
        steps = 20; ds = cap/steps; T = T0; y_acc = np.zeros_like(k)
        for _ in range(steps):
            Rrem = float(csda_range_cm(T))
            T_new = max(T*(1 - ds/max(Rrem, 1e-9)), 1e-3)   # dE/dx ~ T/R scaling
            Tt = T + ME
            y_acc += f_bh(k/Tt) * (k < T) * ds / X0_AL
            T = T_new
        out += n_delta * y_acc * att
    return out

def band(spec_fn, lo, hi, n=200):
    k = np.geomspace(lo, hi, n)
    return np.trapezoid(spec_fn(k), np.log(k))

def gas_conv(k_kev, d_cm=0.3):
    return 1 - np.exp(-loglog(AR_EN, k_kev)*RHO_AR*d_cm)

def run():
    rows = []
    bands = [(5, 30), (30, 150)]
    print(f"{'E[MeV]':>6} {'inc':>4} {'t[mm]':>5} | "
          f"{'prim 5-30':>10} {'prim 30-150':>11} | {'delta 5-30':>10} {'delta 30-150':>12} | "
          f"{'gas-conv/e- (3mm Ar)':>21}")
    for E in (119.0, 144.0):
        for inc in (10, 40):
            for t_mm in (0.5, 1, 2, 3, 5):
                s = t_mm/10/math.cos(math.radians(inc))
                p = [band(lambda k: primary_spectrum(k, E, s), *b) for b in bands]
                d = [band(lambda k: delta_spectrum(k, E, s), *b) for b in bands]
                kk = np.geomspace(5, 150, 200)
                tot = primary_spectrum(kk, E, s) + delta_spectrum(kk, E, s)
                conv = np.trapezoid(tot*gas_conv(kk), np.log(kk))
                rows.append((E, inc, t_mm, *p, *d, conv))
                print(f"{E:6.0f} {inc:4d} {t_mm:5.1f} | {p[0]:10.2e} {p[1]:11.2e} | "
                      f"{d[0]:10.2e} {d[1]:12.2e} | {conv:21.2e}")
    return rows

def plots():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    k = np.geomspace(3, 150, 300)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    cols = ["#1b9e77", "#d95f02", "#7570b3"]
    for c, t in zip(cols, (1, 2, 5)):
        s = t/10/math.cos(math.radians(10))
        ax[0].loglog(k, primary_spectrum(k, 119, s), c=c, label=f"{t} mm primary")
        ax[0].loglog(k, delta_spectrum(k, 119, s), c=c, ls="--", label=f"{t} mm delta-ray")
    ax[0].set(xlabel="photon energy k [keV]", ylabel="photons per electron per ln k",
              title="Photons leaving the Al wall (119 MeV e-, 10° incidence)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which="both")
    ax[1].loglog(k, gas_conv(k), c="k")
    ax[1].set(xlabel="photon energy k [keV]", ylabel="P(interact in 3 mm Ar, 1 atm)",
              title="Gas-gap response to those photons")
    ax[1].grid(alpha=.3, which="both")
    fig.savefig(OUT/"chapeau_photons.png", dpi=150)

if __name__ == "__main__":
    run(); plots()
    # sanity prints
    print("\nX0(Al) = %.2f cm ; MS angle 119 MeV, 2 mm: %.2f deg" % (
        X0_AL, math.degrees(13.6/119*math.sqrt(0.2/X0_AL)*(1+0.038*math.log(0.2/X0_AL)))))
    print("MIP-like dE in 3 mm Ar ~0.76 keV; Ar absorbs 20 keV photon w.p. %.4f, 50 keV %.5f"
          % (gas_conv(20.0), gas_conv(50.0)))
