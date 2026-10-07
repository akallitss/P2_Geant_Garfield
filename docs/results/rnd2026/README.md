# rnd2026 — first look (Geant4 only)

Batch: `submit_condor_p2.py --scan rnd2026`, HTCondor cluster 12403408,
commit f643863, 106 jobs, 80.44 M events, `/eos/user/a/akallits/p2_sim/rnd2026`.
Analysis: `scripts/analysis/rnd2026_photons_electrons.py` (run on lxplus).

Primary-ionization level only: "detected" = ≥ 1 primary electron in the drift
gas; observable = energy deposited in the drift gas. No gain, threshold or
electronics yet (needs the Magboltz tables and Stage B/C).

| result | 3.96 mm drift gas (4 mm frame) |
|---|---|
| P(detect) per photon, 40–100 keV | Ar 2.5 → 0.6 × 10⁻³ · Ne 1.5 → 0.45 × 10⁻³ |
| Ar / Ne | 1.7× at 40 keV, 1.3× at 100 keV — wall conversions (mostly the mesh) dominate |
| isobutane 5 → 10 % | no change in P(detect) or in the e⁻ deposit |
| 5 mm drift gas | P(detect) +5–15 %, e⁻ deposit +40 % |
| e⁻ 100 MeV most probable deposit | Ar 0.8 keV · Ne 0.4 keV |
| detected photons | median 6–10 keV; mesh Fe/Ni K lines (6.4 / 7.5 keV) + full absorption |
| upper cut at 99 % e⁻ eff. | Ar 9–13 keV, Ne 6–8 keV → removes 35–55 % (Ar), 40–70 % (Ne) of detected photons |
| upper cut at 95 % e⁻ eff. | 3–5 keV → removes 70–93 % |
| MM0 rate, P2Sim flux 10–150 keV | Ar 79 → 41 MHz · Ne 37 → 20 MHz after the 99 % cut (whole layer) |

Files: `photon_summary.csv`, `electron_summary.csv`, `fold_p2sim.csv`, `fig1–4`.

Caveats
- Gas vs wall is split on `convVolume`. Up to f643863 `interactionClass` tagged
  a soft (< 20 keV) beam photon absorbed in the gas as fluorescence (fixed in
  `EventAction.cc` after this batch).
- 10–150 keV covers 84 % of the MM0 flux; the softer part is not included yet.
- The 99 % electron cut rests on 10⁴ electrons per point (±~1 keV).
