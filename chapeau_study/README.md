# Chapeau (Al) low-energy-photon study — analytic, no Geant4

Question: can ~119-144 MeV back-scattered electrons (P2 155 MeV beam, elastic
e-p at 140-150 deg -> 119 MeV; 200 MeV option -> 144 MeV; see
`../docs/P2_EXPERIMENT.md` s3) crossing an aluminium chapeau wall make the
low-energy photons seen in the sim?

Run: `python chapeau_photons.py` -> table on stdout + `chapeau_photons.png`.

## Assumptions (nothing about the chapeau is in this repo)
- Al wall, thickness **unknown** (scanned 0.5-5 mm), incidence 10 and 40 deg.
- Sources: primary Bethe-Heitler brems + brems of knock-on electrons. Compton/
  photoabsorption of the hard brems photons, Al fluorescence (1.5 keV, dies in
  microns), TR not modelled. Delta-ray part is good to ~factor 2.
- NIST XCOM attenuation (Al, Ar). Gas gap taken as 3 mm pure Ar, 1 atm.

## Result (per incident electron, photons that exit the wall)
| t (Al) | 5-30 keV | 30-150 keV |
|---|---|---|
| 1 mm | ~1e-2 | ~2.4e-2 |
| 2 mm | ~1.3e-2 | ~4.7e-2 |
| 5 mm | ~1.7e-2 | ~1.1e-1 |
- Scales ~ t/X0 (X0(Al)=8.9 cm); delta-ray brems is a <10% correction.
- Energy-independent between 119 and 144 MeV (brems ~ 1/k, E only sets the cutoff).
- Below ~10 keV the wall itself absorbs them (mu*rho*1 mm ~ 70 at 10 keV).
- Photons are collinear with the electron (MS cone ~1 deg for 2 mm), i.e. they
  arrive with the electron in the same pad/time, not as separate hits.
- Conversion in the 3 mm gas gap: ~1e-4 per electron -> negligible vs. a MIP.

## Caveat / follow-up
Plausible yes as a *photon-count* source (a few % of electrons carry a 30-150 keV
photon), but this is far below the ~1e3 photon/electron design background and the
yield is a smooth 1/k, so it would not make a peak. To compare we need the sim's
actual low-energy photon spectrum + production-process/volume tags, and the real
chapeau thickness. Real Geant4 run (e- at 119 MeV into Al slab + gas, tag by
creator process) is the follow-up.

## Wall thickness estimate (`wall_estimate.py`) and yield plots (`plot_yield.py`)
Freestanding flat Al plate, clamped, Al 6061-T6, SF 3, deflection limit 1-3 mm
(gap is only 3 mm). Span and overpressure are unknown: ~10 mbar over a
150-300 mm half-span gives 0.5-1 mm (stress) to 1-3 mm (stiffness). A real
1 bar-rated wide window would be 5-12+ mm. **Nominal 2 mm, band 1-5 mm.**
Plots: `yield_vs_energy.png` (1/2/5 mm), `yield_by_process.png` (2 mm,
primary vs delta-ray brems, produced vs exiting).
