#!/usr/bin/env python3
"""
Analytic photon-interaction budget for the P2 wedge at 50-150 keV.

Purpose: a back-of-envelope that the Geant4 campaign must reproduce (or
contradict) before anything is trusted.  It answers three questions the
gas-choice argument turns on:

  1. How different are argon and neon in *interaction probability* in the
     drift gap at 50-150 keV?  (Well above the Ar K-edge, so not very.)
  2. How much energy does each interaction channel actually leave in a
     3 mm gap, versus a MIP?  (The discrimination premise.)
  3. How does gas conversion compare with conversions in the *walls* --
     copper pad plane, mesh, mylar -- which inject electrons into the gas
     no matter what the gas is?  (The floor on any gas-based rejection.)

Everything here is deliberately crude: uniform-depth interactions, a
straight-line escape model for photo/Compton electrons, Katz-Penfold
ranges.  Order-of-magnitude only -- the point is to know what the
simulation should come out near, and which channels must be recorded
separately.

Data sources
------------
mu/rho, mu_en/rho : NIST X-Ray Mass Attenuation Coefficients (Hubbell &
    Seltzer), https://physics.nist.gov/PhysRefData/XrayMassCoef/
    element tables z10 (Ne), z18 (Ar), z26 (Fe), z29 (Cu) and the
    compound table for Mylar.  Values transcribed below.
Klein-Nishina : integrated numerically from the differential cross
    section (free-electron approximation -- adequate for Ne/Ar/mylar
    above ~30 keV, optimistic for Cu/Fe where binding matters more).
Electron ranges : Katz & Penfold, Rev. Mod. Phys. 24 (1952) 28,
    R_csda [g/cm2] = 0.412 E^(1.265 - 0.0954 ln E), E in MeV,
    valid 0.01-3 MeV.  Practical penetration depth is taken as
    DETOUR * R_csda.

Usage:  python3 scripts/photon_budget.py
"""

import math

# ----------------------------------------------------------------------
# Constants and knobs
# ----------------------------------------------------------------------

MEC2_KEV = 510.999
R_E_CM = 2.8179403262e-13
N_A = 6.02214076e23
DETOUR = 0.5          # practical penetration / CSDA path length, low-Z, ~50 keV
T_KELVIN = 293.15
P_PA = 101325.0
R_GAS = 8.314462

# NIST mass attenuation coefficients [cm2/g], keyed by energy [keV].
# (mu/rho with coherent, mu_en/rho)
NIST = {
    "Ne": {50: (0.2579, 0.08182), 60: (0.2161, 0.05287),
           80: (0.1781, 0.03273), 100: (0.1600, 0.02733),
           150: (0.1370, 0.02590)},
    "Ar": {50: (0.7012, 0.5020), 60: (0.4664, 0.2904),
           80: (0.2760, 0.1280), 100: (0.2043, 0.07344),
           150: (0.1427, 0.03703)},
    "Cu": {50: (2.613, 2.192), 60: (1.593, 1.290),
           80: (0.7630, 0.5581), 100: (0.4584, 0.2949),
           150: (0.2217, 0.1027)},
    "Mylar": {50: (0.2020, 0.03082), 60: (0.1868, 0.02508),
              80: (0.1695, 0.02247), 100: (0.1586, 0.02297),
              150: (0.1406, 0.02567)},
}

# Molar masses [g/mol] and Z/A for the gases we need.
GAS = {"Ne": dict(M=20.1797, Z=10, A=20.1797),
       "Ar": dict(M=39.948, Z=18, A=39.948)}

# Condensed materials: density [g/cm3], Z, A (effective for scaling).
SOLID = {"Cu": dict(rho=8.96, Z=29, A=63.546),
         "Fe": dict(rho=8.00, Z=26, A=55.845),   # stainless ~ Fe for our purposes
         "Mylar": dict(rho=1.39, Z=None, A=None)}

ENERGIES = [50, 60, 80, 100, 150]


# ----------------------------------------------------------------------
# Klein-Nishina
# ----------------------------------------------------------------------

def klein_nishina(e_kev, n=20000):
    """Return (sigma_tot [cm2/electron], <T> [keV], T_max [keV]).

    Numerical integration of dsigma/dOmega over scattering angle for a
    free electron at rest.
    """
    alpha = e_kev / MEC2_KEV
    sig = 0.0
    sig_t = 0.0
    for i in range(n):
        th = (i + 0.5) * math.pi / n
        dth = math.pi / n
        ratio = 1.0 / (1.0 + alpha * (1.0 - math.cos(th)))
        dsig = 0.5 * R_E_CM**2 * ratio**2 * (ratio + 1.0 / ratio - math.sin(th)**2)
        dom = 2.0 * math.pi * math.sin(th) * dth
        sig += dsig * dom
        sig_t += dsig * dom * e_kev * (1.0 - ratio)
    t_max = e_kev * 2.0 * alpha / (1.0 + 2.0 * alpha)
    return sig, sig_t / sig, t_max


def gas_density(name):
    """Ideal-gas density [g/cm3] at 1 atm, 20 C."""
    return P_PA * GAS[name]["M"] * 1e-3 / (R_GAS * T_KELVIN) * 1e-3


def compton_mu_over_rho(z, a, e_kev):
    """Incoherent mu/rho [cm2/g] in the free-electron approximation."""
    sig, _, _ = klein_nishina(e_kev)
    return sig * N_A * z / a


def csda_range(e_kev):
    """Katz-Penfold CSDA range [g/cm2] for an electron of energy e_kev."""
    e = e_kev / 1000.0
    n = 1.265 - 0.0954 * math.log(e)
    return 0.412 * e**n


def stopping_power(e_kev, de=0.2):
    """Mean collision stopping power [MeV cm2/g] near e_kev, from dR/dE."""
    e1, e2 = e_kev * (1 - de), e_kev * (1 + de)
    return (e2 - e1) / 1000.0 / (csda_range(e2) - csda_range(e1))


# ----------------------------------------------------------------------
# 1. Gas: interaction probability in the drift gap
# ----------------------------------------------------------------------

def gas_table(gap_mm=3.0):
    print(f"\n=== 1. Pure-gas interaction probability, {gap_mm:.1f} mm gap, "
          "1 atm / 20 C ===")
    print("   (mu/rho from NIST; photoelectric = total - incoherent - "
          "coherent, with")
    print("    incoherent from Klein-Nishina.  Negative/zero PE at high E "
          "means the")
    print("    coherent part is being absorbed into the residual -- treat "
          "as '~0'.)\n")
    hdr = (f"{'E[keV]':>7} | {'P_int(Ar)':>10} {'P_int(Ne)':>10} "
           f"{'ratio':>6} | {'f_PE(Ar)':>9} {'f_PE(Ne)':>9} | "
           f"{'<E_tr|int> Ar':>13} {'Ne':>7}")
    print(hdr)
    print("-" * len(hdr))
    out = {}
    for e in ENERGIES:
        row = {}
        for g in ("Ar", "Ne"):
            mu_r, muen_r = NIST[g][e]
            rho = gas_density(g)
            p = 1.0 - math.exp(-mu_r * rho * gap_mm / 10.0)
            comp = compton_mu_over_rho(GAS[g]["Z"], GAS[g]["A"], e)
            f_pe = max(0.0, (mu_r - comp)) / mu_r
            row[g] = dict(p=p, f_pe=f_pe, e_tr=muen_r / mu_r * e)
        out[e] = row
        print(f"{e:>7} | {row['Ar']['p']:>10.3e} {row['Ne']['p']:>10.3e} "
              f"{row['Ar']['p']/row['Ne']['p']:>6.2f} | "
              f"{row['Ar']['f_pe']:>8.1%} {row['Ne']['f_pe']:>9.1%} | "
              f"{row['Ar']['e_tr']:>10.1f} keV {row['Ne']['e_tr']:>4.1f} keV")
    print("\n  f_PE here lumps coherent scattering in with 'not Compton', so "
          "it is an\n  upper bound on the true photoelectric fraction "
          "(coherent deposits nothing).")
    return out


# ----------------------------------------------------------------------
# 2. What each channel deposits in the gap, versus a MIP
# ----------------------------------------------------------------------

MIP_DEDX = {"Ar": 2.53, "Ne": 1.45}      # keV/cm, PDG 2024 Table 35.5, pure gas


def deposit_bracket(t_kev, rho, gap_mm):
    """Bracket the energy an electron of energy t_kev leaves in the gap.

    Lower bound: straight-line crossing at the entry-energy stopping power
    (right for range >> gap, ignores the rise of dE/dx as it slows).
    Upper bound: min(T, T * gap / R_p) -- fraction-of-practical-range
    (right near containment, over-counts for a fast crossing electron).
    Returns (lo, hi, R_p [mm]).
    """
    rp = DETOUR * csda_range(t_kev) / rho * 10.0          # mm
    lo = min(t_kev, stopping_power(t_kev) * 1000.0 * rho * gap_mm / 10.0)
    hi = t_kev * min(1.0, gap_mm / rp)
    return min(lo, hi), max(lo, hi), rp


def deposit_table(gap_mm=3.0):
    print(f"\n=== 2. Energy left in a {gap_mm:.1f} mm gap, per channel ===")
    print("    Deposits are bracketed [dE/dx crossing .. fraction-of-range];")
    print("    the two bounds converge when the electron is contained.")
    for g in ("Ar", "Ne"):
        rho = gas_density(g)
        mip = MIP_DEDX[g] * gap_mm / 10.0
        print(f"\n  --- {g} (rho = {rho*1e3:.3f} mg/cm3, "
              f"MIP mean edep = {mip*1000:.0f} eV over {gap_mm:.0f} mm) ---")
        hdr = (f"  {'E[keV]':>7} | {'T_max':>6} {'<T>':>6} {'R_p':>7} "
               f"{'Compton dep [keV]':>18} {'xMIP':>10} | "
               f"{'R_p(PE)':>8} {'PE dep [keV]':>14} {'xMIP':>10}")
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        for e in ENERGIES:
            _, t_mean, t_max = klein_nishina(e)
            c_lo, c_hi, rp_c = deposit_bracket(t_mean, rho, gap_mm)
            # Photoelectron carries ~the full photon energy (Ar K-fluorescence
            # escape removes 2.96 keV in 12 % of Ar photoabsorptions).
            p_lo, p_hi, rp_p = deposit_bracket(e, rho, gap_mm)
            print(f"  {e:>7} | {t_max:>5.1f}k {t_mean:>5.1f}k {rp_c:>6.2f}mm "
                  f"{c_lo:>8.2f}-{c_hi:<9.2f} {c_lo/mip:>4.1f}-{c_hi/mip:<5.1f} | "
                  f"{rp_p:>6.1f}mm {p_lo:>6.2f}-{p_hi:<7.2f} "
                  f"{p_lo/mip:>4.1f}-{p_hi/mip:<5.1f}")
        print("    Means only.  The Compton spectrum is broad and reaches 0, "
              "so a tail of\n    Compton events lands at or below the MIP "
              "deposit -- that tail, plus the\n    wall-crossing electrons of "
              "table 3, is what limits a charge cut.")


# ----------------------------------------------------------------------
# 3. Wall conversions -- the floor no gas choice can move
# ----------------------------------------------------------------------

def escape_fraction_slab(thick_um, range_um):
    """Fraction of uniformly-created electrons escaping one face of a slab."""
    return 0.5 * min(1.0, range_um / thick_um)


def escape_fraction_wire(diam_um, range_um):
    """Same for a cylindrical wire: shell within one range of the surface."""
    r = diam_um / 2.0
    if range_um >= r:
        return 0.5
    return 0.5 * (1.0 - ((r - range_um) / r)**2)


def wall_table(gap_mm=3.0, e_kev=60):
    print(f"\n=== 3. Electron injection into the gas at {e_kev} keV "
          f"(per incident photon) ===")
    print("    Layers as built in DetectorConstruction::ConstructP2().\n")

    # Compton and photoelectric split per material at this energy.
    def split(matname, z, a):
        mu_r, _ = NIST[matname][e_kev]
        comp = compton_mu_over_rho(z, a, e_kev) if z else None
        if comp is None:                     # mylar: Z/A ~ 0.520
            comp = compton_mu_over_rho(0.520, 1.0, e_kev)
        f_c = min(1.0, comp / mu_r)
        return mu_r, f_c, 1.0 - f_c

    _, t_mean, _ = klein_nishina(e_kev)
    rows = []

    # -- Copper pad plane (18 um, 98.3 % coverage), faces the amp gap.
    mu_r, f_c, f_pe = split("Cu", 29, 63.546)
    rho = SOLID["Cu"]["rho"]
    t_um = 18.0 * 0.983
    p_int = 1.0 - math.exp(-mu_r * rho * t_um * 1e-4)
    rp_pe = DETOUR * csda_range(e_kev) / rho * 1e4          # um
    rp_c = DETOUR * csda_range(t_mean) / rho * 1e4
    inj = p_int * (f_pe * escape_fraction_slab(t_um, rp_pe)
                   + f_c * escape_fraction_slab(t_um, rp_c))
    rows.append(("PCB_Cu_F (18 um Cu, -> amp gap)", p_int, f_pe, rp_pe, inj))

    # -- Micromesh: 19 um stainless wires, 48 um opening.  Fe scaled from Cu.
    mu_cu, f_c_cu, f_pe_cu = split("Cu", 29, 63.546)
    comp_fe = compton_mu_over_rho(26, 55.845, e_kev)
    pe_fe = (mu_cu - compton_mu_over_rho(29, 63.546, e_kev)) \
        * (26 / 29)**4.5 * (63.546 / 55.845)
    mu_fe = comp_fe + pe_fe
    f_pe_fe = pe_fe / mu_fe
    rho = SOLID["Fe"]["rho"]
    # Solid-equivalent areal thickness of a plain-weave square mesh: one x-
    # and one y-wire per pitch^2 cell, so t_eff = 2*(pi d^2/4)/pitch.
    # (This is exactly what DetectorConstruction builds: a 2*d slab at fill
    # fraction pi*d/(4*pitch) has the same areal density -- verified. An
    # earlier version of this script used a projected-coverage estimate that
    # was 1.7x too high, which inflated the mesh conversion rate.)
    d_wire, pitch = 19.0, 67.0
    t_eff_um = 2.0 * (math.pi * d_wire**2 / 4.0) / pitch
    p_int = 1.0 - math.exp(-mu_fe * rho * t_eff_um * 1e-4)
    rp_pe = DETOUR * csda_range(e_kev) / rho * 1e4
    inj = p_int * (f_pe_fe * escape_fraction_wire(19.0, rp_pe)
                   + (1 - f_pe_fe) * escape_fraction_wire(
                       19.0, DETOUR * csda_range(t_mean) / rho * 1e4))
    rows.append(("Micromesh (19 um SS wire, -> both gaps)",
                 p_int, f_pe_fe, rp_pe, inj))

    # -- Drift cathode: 2 x 12 um mylar + 0.1 um Al, and the 40 um window.
    mu_r, f_c, f_pe = split("Mylar", None, None)
    rho = SOLID["Mylar"]["rho"]
    for label, t_um in (("Drift cathode mylar (2 x 12 um)", 24.0),
                        ("Front window mylar (40 um)", 40.0)):
        p_int = 1.0 - math.exp(-mu_r * rho * t_um * 1e-4)
        rp_pe = DETOUR * csda_range(e_kev) / rho * 1e4
        rp_c = DETOUR * csda_range(t_mean) / rho * 1e4
        inj = p_int * (f_pe * escape_fraction_slab(t_um, rp_pe)
                       + f_c * escape_fraction_slab(t_um, rp_c))
        rows.append((label, p_int, f_pe, rp_pe, inj))

    hdr = (f"  {'layer':<40} {'P_int':>9} {'f_PE':>6} {'R_p[um]':>8} "
           f"{'P(e- into gas)':>15}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for name, p, fpe, rp, inj in rows:
        print(f"  {name:<40} {p:>9.2e} {fpe:>6.1%} {rp:>8.1f} {inj:>15.2e}")

    # Gas conversions for comparison.
    print()
    for g in ("Ar", "Ne"):
        mu_r, _ = NIST[g][e_kev]
        p = 1.0 - math.exp(-mu_r * gas_density(g) * gap_mm / 10.0)
        print(f"  {'drift gas, ' + g + f' ({gap_mm:.0f} mm)':<40} "
              f"{p:>9.2e} {'':>6} {'':>8} {p:>15.2e}")

    tot_wall = sum(r[4] for r in rows)
    print(f"\n  wall total ~ {tot_wall:.2e} electrons into gas per photon.")
    for g in ("Ar", "Ne"):
        mu_r, _ = NIST[g][e_kev]
        p = 1.0 - math.exp(-mu_r * gas_density(g) * gap_mm / 10.0)
        print(f"    wall/gas ratio in {g}: {tot_wall/p:>5.1f}")

    # Cu K-fluorescence re-absorbed in the gas: a gas-sensitive secondary.
    print("\n  Cu-K fluorescence channel (8.05 keV photon out of the pad "
          "plane,\n  re-absorbed in the gas -> point-like ~8 keV blob):")
    mu_cu, _ = NIST["Cu"][e_kev]
    p_cu = 1.0 - math.exp(-mu_cu * SOLID["Cu"]["rho"] * 18.0 * 0.983 * 1e-4)
    omega_k = 0.44                       # Cu K fluorescence yield
    f_esc = 0.30                         # 8 keV out of 18 um Cu, rough
    # 8 keV absorption in the gas: extrapolate NIST 10 keV point as E^-3 (PE).
    for g in ("Ar", "Ne"):
        mu8 = {"Ar": 63.16, "Ne": 11.97}[g] * (10.0 / 8.05)**3
        p_abs = 1.0 - math.exp(-mu8 * gas_density(g) * gap_mm / 10.0)
        chain = p_cu * f_pe_cu * omega_k * f_esc * p_abs
        print(f"    {g}: P(8 keV absorbed in {gap_mm:.0f} mm) = {p_abs:.3f}"
              f"  ->  chain = {chain:.2e} per incident photon")


if __name__ == "__main__":
    print(__doc__.split("Usage:")[0].strip()[:0] or "", end="")
    print("P2 photon-interaction budget -- analytic cross-check for the "
          "Geant4 campaign")
    print("=" * 72)
    gas_table()
    deposit_table()
    for e in (60, 100):
        wall_table(e_kev=e)
    print("\nDone.  These are order-of-magnitude targets for the simulation, "
          "not results.")
