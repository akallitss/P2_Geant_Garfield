#!/usr/bin/env python3
"""Rough minimum thickness of a freestanding flat Al wall holding gas.
Clamped circular plate, uniform dp: sigma_max = 0.75 dp (a/t)^2 ; w0 = dp a^4/(64 D),
D = E t^3 / (12 (1-nu^2)).  Al 6061-T6: E=69 GPa, nu=0.33, yield ~240 MPa, SF 3.
Radius a = half-span of the flat window (UNKNOWN; scanned)."""
import numpy as np
E, nu, sy, SF = 69e9, 0.33, 240e6, 3.0
print(f"{'dp':>8} {'a[mm]':>6} | {'t_stress[mm]':>12} {'t for w<1mm':>12} {'t for w<3mm':>12}")
for dp_mbar in (10, 100, 1000):
    dp = dp_mbar*100.0
    for a_mm in (150, 300, 650):
        a = a_mm/1e3
        t_s = a*np.sqrt(0.75*dp*SF/sy)
        def t_w(w): return (dp*a**4*12*(1-nu**2)/(64*E*w))**(1/3)
        print(f"{dp_mbar:6d}mb {a_mm:6d} | {t_s*1e3:12.2f} {t_w(1e-3)*1e3:12.2f} {t_w(3e-3)*1e3:12.2f}")
