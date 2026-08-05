"""Python port of the ATLAS Athena VMM shaper response.

Faithful reimplementation of ``reference/VMM_Shaper.cxx`` (Athena
MM_Digitization, Apache-2.0, constants by G. Iakovidis — described in
CERN-THESIS-2014-148 §7.1.3). The VMM front-end shaper is third-order
"Delayed Dissipative Feedback": one real pole and a complex-conjugate pole
pair. The closed-form response to a delta charge Q arriving at t0 is

    h(t) = Q * S * a^3 * p0 * |p1|^2 *
           [ K0 * exp(-p0 (t-t0))
             + 2|K1| * exp(-Re(p1)(t-t0)) * cos(-Im(p1)(t-t0) + arg K1) ]

with a = t_peak/1.5, p0 = 1.263/a, p1 = (1.149 - 0.786 i)/a,
K0 = 1.584, K1 = -0.792 - 0.115 i, and S chosen such that the peak of h for
a single electron equals its (ballistic-scaled) charge:
S = chargeScaleFactor * peakTimeChargeScaling, chargeScaleFactor = 1/0.411819.

Ballistic deficit: the Micromegas ion tail lasts ~150 ns; for peaking times
below that only a fraction t_peak/150 ns of the charge is integrated
(peakTimeChargeScaling). At t_peak = 200 ns the full charge is seen.

Units: time in ns, charge in whatever unit you feed in (Athena uses
electrons); the returned amplitude is in the same charge unit. Thresholds
are therefore also in charge units.

Differences from the C++ (documented, intentional):
- peak/threshold searches run on a dense numpy time grid instead of
  Athena's adaptive coarse/fine scan — same semantics (first local maximum
  above threshold; last below-threshold sample before it), simpler code.
- everything is vectorized over the time grid.

Validated by ``python vmm_shaper.py`` (self-test, see bottom).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# Constants from reference/VMM_Shaper.cxx (G. Iakovidis)
K0 = 1.584
RE_K1 = -0.792
IM_K1 = -0.115
CHARGE_SCALE_FACTOR = 1.0 / 0.411819
MM_ION_FLOW_TIME_NS = 150.0

ELECTRONS_PER_FC = 6241.509  # 1 fC in electrons


@dataclass
class ShaperOutput:
    """Result of digitizing one channel."""

    fired: bool
    amplitude: float = 0.0   # shaped peak amplitude, input charge units
    peak_time: float = 0.0   # ns, time of first peak above threshold
    threshold_time: float = 0.0  # ns, threshold-crossing time (if computed)


class VMMShaper:
    """VMM shaper response for one channel.

    Parameters
    ----------
    peak_time : shaper peaking time in ns (VMM3a: 25, 50, 100 or 200).
    lower_time_window, upper_time_window : ns, acceptance window for the
        peak / threshold crossing. Athena extends the search window by the
        VMM dead time below and a "graze" window above; do that in the
        caller (see vmm_emulator.PadDigitizer).
    time_step : ns, grid resolution of the search (Athena fine step: 0.1).
    """

    def __init__(self, peak_time: float, lower_time_window: float,
                 upper_time_window: float, time_step: float = 0.1):
        self.peak_time = float(peak_time)
        self.lower_time_window = float(lower_time_window)
        self.upper_time_window = float(upper_time_window)
        self.time_step = float(time_step)

        a = self.peak_time / 1.5  # ns
        self.a = a
        self.pole0 = 1.263 / a
        self.re_pole1 = 1.149 / a
        self.im_pole1 = -0.786 / a
        pole1_square = self.re_pole1 ** 2 + self.im_pole1 ** 2
        self.k1_abs = math.hypot(RE_K1, IM_K1)
        self.arg_k1 = math.atan2(IM_K1, RE_K1)

        # ballistic deficit for peaking times below the ion-tail duration
        self.peak_time_charge_scaling = (
            self.peak_time / MM_ION_FLOW_TIME_NS
            if self.peak_time < MM_ION_FLOW_TIME_NS else 1.0
        )
        # a^3 * p0 * |p1|^2 is dimensionless when poles carry 1/a
        self.pre_calc = (CHARGE_SCALE_FACTOR * self.peak_time_charge_scaling
                         * a ** 3 * self.pole0 * pole1_square)

    # ------------------------------------------------------------------
    def response(self, charges, times, t):
        """Shaped response at time(s) ``t`` (ns) for electrons of
        ``charges`` arriving at ``times``. Vectorized over ``t``."""
        t = np.atleast_1d(np.asarray(t, dtype=float))
        charges = np.asarray(charges, dtype=float)
        times = np.asarray(times, dtype=float)
        dt = t[:, None] - times[None, :]          # (nt, ne)
        active = dt >= 0.0
        dtc = np.where(active, dt, 0.0)
        shape = (K0 * np.exp(-dtc * self.pole0)
                 + 2.0 * self.k1_abs * np.exp(-dtc * self.re_pole1)
                 * np.cos(-dtc * self.im_pole1 + self.arg_k1))
        out = (charges[None, :] * self.pre_calc * np.where(active, shape, 0.0)
               ).sum(axis=1)
        return out if out.size > 1 else float(out[0])

    # ------------------------------------------------------------------
    def _grid(self, times):
        start = max(self.lower_time_window, float(np.min(times)))
        stop = self.upper_time_window
        if stop <= start:
            return None
        n = int(round((stop - start) / self.time_step)) + 1
        return np.linspace(start, stop, n)

    def above_threshold_simple(self, charges, times, threshold) -> bool:
        """Cheap pre-filter: ballistic-scaled charge sum vs threshold
        (Athena aboveThresholdSimple)."""
        charges = np.asarray(charges, dtype=float)
        times = np.asarray(times, dtype=float)
        sel = ((times >= self.lower_time_window - self.peak_time)
               & (times <= self.upper_time_window))
        return charges[sel].sum() * self.peak_time_charge_scaling >= threshold

    def has_charge_above_threshold(self, charges, times, threshold) -> bool:
        if len(charges) == 0 or not self.above_threshold_simple(
                charges, times, threshold):
            return False
        grid = self._grid(times)
        if grid is None:
            return False
        return bool(np.any(self.response(charges, times, grid) >= threshold))

    def peak_response(self, charges, times, threshold) -> ShaperOutput:
        """First local maximum of the shaped signal above threshold
        (Athena vmmPeakResponse). Amplitude in input charge units."""
        if len(charges) == 0:
            return ShaperOutput(False)
        grid = self._grid(times)
        if grid is None:
            return ShaperOutput(False)
        r = self.response(charges, times, grid)
        # local maxima: r[i-1] <= r[i] > r[i+1]
        peaks = np.flatnonzero((r[1:-1] >= r[:-2]) & (r[1:-1] > r[2:])) + 1
        peaks = peaks[r[peaks] >= threshold]
        if peaks.size == 0:
            return ShaperOutput(False)
        i = int(peaks[0])
        return ShaperOutput(True, float(r[i]), float(grid[i]))

    def threshold_response(self, charges, times, threshold) -> ShaperOutput:
        """Peak amplitude + threshold-crossing time
        (Athena vmmThresholdResponse)."""
        out = self.peak_response(charges, times, threshold)
        if not out.fired:
            return out
        grid = self._grid(times)
        r = self.response(charges, times, grid)
        above = np.flatnonzero(r >= threshold)
        out.threshold_time = float(grid[above[0]]) if above.size else out.peak_time
        return out


# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Self-test: normalization, ballistic scaling, peak position.
    q = 10000.0  # electrons
    for tp in (25.0, 50.0, 100.0, 200.0):
        sh = VMMShaper(tp, lower_time_window=-50.0, upper_time_window=600.0)
        out = sh.peak_response([q], [0.0], threshold=0.0)
        expected = q * sh.peak_time_charge_scaling
        assert out.fired
        assert abs(out.amplitude - expected) / expected < 0.01, \
            (tp, out.amplitude, expected)
        assert abs(out.peak_time - tp) < 0.05 * tp + 1.0, (tp, out.peak_time)
    # two electrons, spread in time: amplitude below sum, above single
    sh = VMMShaper(200.0, -50.0, 600.0)
    out2 = sh.threshold_response([q, q], [0.0, 100.0], threshold=q / 2)
    assert q < out2.amplitude < 2 * q
    assert out2.threshold_time < out2.peak_time
    print("vmm_shaper self-test OK")
    print(f"  peak(200ns, 1e4 e-): amp={out2.amplitude:.0f} e-, "
          f"t_peak={out2.peak_time:.1f} ns, t_thr={out2.threshold_time:.1f} ns")
