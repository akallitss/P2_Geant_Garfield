"""gas.py — drift transport parameters, with Magboltz as a drop-in.

Stage B needs, per gas and drift field: drift velocity v_d, longitudinal and
transverse diffusion, and the attachment length. None of those are measured
or computed yet — P0.10 (Magboltz tables) is open — so this module ships
**placeholders** and is built so that swapping them for real tables changes
one function and nothing else. That separation is the point of P0.9's "write
it with pluggable gas tables; unit-test with fake tables so Magboltz becomes
a drop-in".

THREE WARNINGS ABOUT THE NUMBERS BELOW, in decreasing order of how much they
can move a result:

1. **v_d and the diffusion constants are placeholders**, taken from the same
   table as `vmm/time_resolution.py` so the two stay consistent. The drift
   times are the ones quoted in the campaign docs (Ar-based 60-75 ns, Ne/CO2
   ~110 ns over 3 mm); sigma_L is a generic 230-300 um/sqrt(cm). Anything
   that scales as 1/v_d — which is most of the timing — inherits their error
   directly.

2. **Water contamination dominates v_d in a real gas system.** MX17's bench
   measures 36.6 um/ns against a much higher dry-Magboltz prediction in
   Ar/iso and infers 1-2 % H2O at the SPS. So even the real dry table will
   sit above the real detector. P0.10 now generates wet variants; until then,
   `v_scale` exists so a run can state which assumption it was made under,
   and it is recorded in the output. Do NOT bake a preference in here — that
   is MX17's explicit advice after doing it the other way.

3. **Attachment is a placeholder of "none"** for the argon mixtures and a
   token value for the CO2-bearing ones. Over 3 mm it is a sub-percent effect
   at any sane oxygen level, which is why it is safe to leave crude — but it
   is not zero-risk if the gas system leaks.

The `TransportTable` interface is what a Magboltz loader must satisfy. When
P0.10 lands, add `MagboltzTable(path)` implementing the same three methods
and pass it to the digitizer; nothing else changes.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np


# ----------------------------------------------------------------------
# Placeholder transport parameters, keyed by the campaign gas names used in
# include/GasMixtures.hh. Values shared with vmm/time_resolution.py.
#
#   n_p       primary clusters / cm          (PDG 2024 Table 35.5, vol-weighted)
#   n_tot     total ionization electrons / cm
#   t_3mm     full-gap drift time at 3 mm [ns]  -> v_d = 3 / t_3mm
#   sigma_L   longitudinal diffusion [um/sqrt(cm)]
#   sigma_T   transverse diffusion   [um/sqrt(cm)]
#   lambda_a  attachment length [cm]; inf = no attachment
# ----------------------------------------------------------------------
_PLACEHOLDER = {
    "ArIso":         dict(n_p=28.3, n_tot=103.2, t_3mm=60.0,  sigma_L=250.0, sigma_T=270.0, lambda_a=np.inf),
    "ArCO2Iso9352":  dict(n_p=26.8, n_tot=99.6,  t_3mm=65.0,  sigma_L=230.0, sigma_T=250.0, lambda_a=500.0),
    "NeIso9010":     dict(n_p=20.7, n_tot=58.0,  t_3mm=75.0,  sigma_L=250.0, sigma_T=280.0, lambda_a=np.inf),
    "NeIso8020":     dict(n_p=24.0, n_tot=70.0,  t_3mm=70.0,  sigma_L=240.0, sigma_T=265.0, lambda_a=np.inf),
    "NeCO2Iso9532":  dict(n_p=15.2, n_tot=45.4,  t_3mm=110.0, sigma_L=300.0, sigma_T=320.0, lambda_a=400.0),
    "NeIsoCO29055":  dict(n_p=18.5, n_tot=53.0,  t_3mm=95.0,  sigma_L=280.0, sigma_T=300.0, lambda_a=450.0),
    "NeCH4937":      dict(n_p=14.1, n_tot=41.0,  t_3mm=70.0,  sigma_L=280.0, sigma_T=300.0, lambda_a=np.inf),
    "NeC2H69010":    dict(n_p=17.0, n_tot=50.0,  t_3mm=72.0,  sigma_L=270.0, sigma_T=290.0, lambda_a=np.inf),
}

# Inherited MX17 mixtures fall back to argon-like transport so that a run in
# one of them does not crash; flagged in `provenance` so it cannot be quoted
# by accident.
_FALLBACK = _PLACEHOLDER["ArIso"]


@dataclass
class Transport:
    """Transport parameters for one (gas, field) point."""

    gas: str
    v_drift_mm_per_ns: float
    sigma_L_mm_per_sqrt_mm: float
    sigma_T_mm_per_sqrt_mm: float
    attach_length_mm: float
    n_p_per_cm: float
    n_tot_per_cm: float
    v_scale: float
    provenance: str

    def to_dict(self):
        d = asdict(self)
        # inf does not survive JSON
        if not np.isfinite(d["attach_length_mm"]):
            d["attach_length_mm"] = -1.0
        return d


class TransportTable:
    """Interface a gas-parameter source must satisfy.

    A Magboltz loader (P0.10) implements this and is passed to the digitizer
    in place of PlaceholderTable. Keep the signature: `for_gas` is the only
    thing Stage B calls.
    """

    def for_gas(self, gas: str, drift_field_V_per_cm: float | None = None) -> Transport:
        raise NotImplementedError


class PlaceholderTable(TransportTable):
    """The table above. Field-independent, which is itself a placeholder.

    v_scale multiplies the drift velocity so a run can be made under an
    explicit wet-gas assumption (e.g. 0.94 for the MX17 bench ratio) and say
    so in its output, rather than a preference being hidden in a constant.
    """

    def __init__(self, v_scale: float = 1.0):
        self.v_scale = float(v_scale)

    def for_gas(self, gas: str, drift_field_V_per_cm: float | None = None) -> Transport:
        known = gas in _PLACEHOLDER
        p = _PLACEHOLDER.get(gas, _FALLBACK)
        prov = ("PLACEHOLDER (campaign doc drift times, PDG n_p); "
                "NOT Magboltz — P0.10 open")
        if not known:
            prov = (f"PLACEHOLDER, and '{gas}' is NOT a campaign gas: argon-like "
                    f"transport substituted. Do not quote timing from this run.")
        v = (3.0 / p["t_3mm"]) * self.v_scale          # mm/ns over a 3 mm gap
        return Transport(
            gas=gas,
            v_drift_mm_per_ns=v,
            # table is um/sqrt(cm); we want mm/sqrt(mm):
            #   sigma[mm] = c[um/sqrt(cm)] * 1e-3 * sqrt(z[mm]/10)
            sigma_L_mm_per_sqrt_mm=p["sigma_L"] * 1e-3 / np.sqrt(10.0),
            sigma_T_mm_per_sqrt_mm=p["sigma_T"] * 1e-3 / np.sqrt(10.0),
            attach_length_mm=(p["lambda_a"] * 10.0
                              if np.isfinite(p["lambda_a"]) else np.inf),
            n_p_per_cm=p["n_p"],
            n_tot_per_cm=p["n_tot"],
            v_scale=self.v_scale,
            provenance=prov,
        )


class JSONTable(TransportTable):
    """Load transport parameters from a JSON file.

    The escape hatch for trying a table without touching code — and the shape
    a Magboltz export should be written into. One object per gas:

        {"NeIso9010": {"v_drift_mm_per_ns": 0.041,
                       "sigma_L_um_sqrtcm": 250, "sigma_T_um_sqrtcm": 280,
                       "attach_length_cm": null,
                       "n_p_per_cm": 20.7, "n_tot_per_cm": 58.0}}
    """

    def __init__(self, path, v_scale: float = 1.0):
        self.path = Path(path)
        self.v_scale = float(v_scale)
        self._d = json.loads(self.path.read_text())

    def for_gas(self, gas: str, drift_field_V_per_cm: float | None = None) -> Transport:
        if gas not in self._d:
            raise KeyError(f"{self.path} has no entry for gas '{gas}' "
                           f"(has: {sorted(self._d)})")
        p = self._d[gas]
        la = p.get("attach_length_cm")
        return Transport(
            gas=gas,
            v_drift_mm_per_ns=float(p["v_drift_mm_per_ns"]) * self.v_scale,
            sigma_L_mm_per_sqrt_mm=float(p["sigma_L_um_sqrtcm"]) * 1e-3 / np.sqrt(10.0),
            sigma_T_mm_per_sqrt_mm=float(p["sigma_T_um_sqrtcm"]) * 1e-3 / np.sqrt(10.0),
            attach_length_mm=(np.inf if la in (None, 0) else float(la) * 10.0),
            n_p_per_cm=float(p.get("n_p_per_cm", np.nan)),
            n_tot_per_cm=float(p.get("n_tot_per_cm", np.nan)),
            v_scale=self.v_scale,
            provenance=f"JSON table {self.path.name}",
        )
