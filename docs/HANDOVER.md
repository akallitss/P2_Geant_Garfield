# Handover: P2 Micromegas simulation

Written 2026-11-01 by Alexandra Kallitsopoulou (CEA Saclay) for the student taking over.
State as of 2026-10-21. Read this first, then `SETUP.md` to get running and
`../CLAUDE.md` for the working rules. `HOWTO.md` has recipes for everyday tasks once you are
set up, including setting up your own accounts and environment from nothing (§0).

## 1. The question

P2 at MESA (Mainz) uses wedge-shaped Micromegas with 1280 pads read out by VMM3a. At the
July 2026 SPS test the VMM was **85 %** efficient, against **96 % for DREAM** on the same pads.
The gain varies by about ×4 across the chamber, so low-gain pads fall under the VMM threshold.

**Goal:** find how to **raise the gas gain** to make up for it (isobutane 5 → 10 %, drift gap
4 → 5 mm, or Ne/iso), while keeping **time resolution ≲ 20 ns** and, if possible, e/γ
separation. The VMM threshold itself is being lowered separately in lab tests.

The end question: **how do drift gap and gas change efficiency and timing resolution, and
which of the two matters more for track reconstruction, in numbers?** The simulation is
reliable on the drift side (primaries, transparency, timing) and much weaker on the
amplification-gap side (gain), which must be calibrated to data. The full set of questions,
and what can answer each, is in `OPEN_PROBLEMS.md` §0.

Particles:
- signal: 100 MeV electrons
- background: 40–100 keV photons, from Matthieu Gailliard's background simulation
- SPS: 150 GeV muons, used only to compare with data

## 2. The chain and where it stands

| step | code | state |
|---|---|---|
| ionization in the gas | `mm_sim` (Geant4 + Heed) | done and validated. Production sample `rnd2026_heed` (106 points) |
| mesh transparency | `scripts/garfield/` + COMSOL maps | done, validated against published data (cold gas only) |
| gain | `scripts/garfield/mesh_gain.cc` | Penning r fitted to data (≈ 0.36) |
| drift → pads | `stage_b/` | runs, **but transparency and pillar losses are not in it yet** |
| electronics | `vmm/` | runs standalone, **not yet on the Stage B output** |

The drift-gap and gas question doesn't need Stage B. The timing toy `vmm/time_resolution.py`
already goes from primaries to VMM efficiency and σ_t, with the gap as a parameter. It needs
updating with the validated inputs (Magboltz, Heed, Garfield transparency), and that update
is the main task (`OPEN_PROBLEMS.md` §2, item 5a). Stage B → VMM is the later, complete
check. One discrepancy with data also needs understanding (§4).

## 3. Numbers to know

Mesh transparency at the SPS point (mesh −450 V, drift 625 V/cm, measured mesh 18 µm wire /
50 µm opening): Ar/iso 95/5 **94.3 %**, Ar/iso 90/10 99.2 %, Ne/iso 95/5 93.8 %,
Ne/iso 85/15 99.9 %, Ar/CO2/iso 93/5/2 97.5 %, Ar/CF4/iso 88/10/2 95.5 %
(`results/mesh/README.md`). These are single runs of a few thousand electrons, so each is
±0.4 % or worse. A 20 000-electron rerun of Ar/iso 95/5 gives 93.4 ± 0.2 % (`OPEN_PROBLEMS.md` §2,
item 7). The losses are electrons diffusing onto the wires, so mixtures with
"hotter" electrons lose more.

Production Geant4 physics is `--physics heed` (Heed for the gas, no PIXE). Older samples
(`rnd2026`, `rnd2026_mylar50`) overcount the primaries. **Don't use them for efficiency.**

**SPS data trap:** take the gas from the corrected `run_config.json` only:
- runs ≤ 58: **Ar/CO2/iso 93/5/2**
- runs 59–72: **Ar/CF4/iso 88/10/2**

**There was no Ar/iso 95/5 at SPS**; that's the cosmic bench. The beam logbook says otherwise
and is wrong. DREAM runs were at drift/mesh 700/450 V, VMM runs at 750/450 V.

## 4. The open problem in one paragraph

Simulated transparency starts to drop ~150–200 V/cm earlier than the data do, and the SPS
data also lose signal at low drift field. An attachment (contamination) model fits the whole
SPS scan, but it needs an implausible amount of O2. A set of geometry variations (gap, wire,
spacing, voltage) has been simulated but **not yet compared to data**. Details, numbers and
the full to-do list: `OPEN_PROBLEMS.md`.

## 5. Your first weeks

1. **Week 1: set up and reproduce.** Follow `SETUP.md`. Run its checks, which reproduce
   a published transparency point and a production Geant4 point. Don't change anything until
   they pass. Fix any step in the docs that was wrong.
2. **Read:** this file, `../CLAUDE.md`, `results/mesh/README.md`, then the logbook slides.
3. **First task: geometry vs data** (`OPEN_PROBLEMS.md` §2, item 1). Take the finished
   geometry runs (`$P2_DATA/condor_transp/`), plot each variation's transparency vs drift field
   against the SPS and cosmic drift scans, and answer: *can a plausible geometry change explain
   the 150–200 V/cm shift?* It needs no new simulation, so it's a good way into the data, the
   scripts and the rules.
4. Then **the drift-gap and gas answer** (`OPEN_PROBLEMS.md` §2, item 5a): update the timing
   toy with the validated inputs, check it against the SPS timing, and compare gaps and gases
   at equal gain. This is the deliverable.

## 6. Where things are

- **Code and docs:** this repo. `../README.md` has the map.
- **Data** (field maps, gas tables, production samples): the data bundle (`SETUP.md` §2).
  Paths written `/eos/user/a/akallits/p2_sim/...` are `$P2_DATA/...` for you.
- **SPS and cosmic data analysis:** the separate `P2_basket_analysis` repo.
- **Logbook** (one dated slide per step):
  https://claude.ai/artifact/TmVumv928tiDApPxyDqSaA (ask Alexandra to share it). Its source is in
  `logbook/project/`.
- **Design files:** `../design/` (gerbers, masks). Every geometry number and its source:
  `P2_GEOMETRY.md`. What the simulation builds: `P2_MODEL.md`.

## 7. People

- **Alexandra Kallitsopoulou** (CEA Saclay): geometry, run conditions, history.
- **Matthieu Gailliard**: P2Sim and the background photon spectra (`BkgHistograms.root`).
- **Dylan Neff**: MX17 detectors. He uses this model and owns the SPS run_config files.
- The Saclay colleague who built the COMSOL woven-mesh model: needed for any new field solve.
