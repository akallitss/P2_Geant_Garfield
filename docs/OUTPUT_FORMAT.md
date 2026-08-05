# Simulation output format (p2 mode)

**Written:** 2026-08-05, after the provenance/scoring rework. Supersedes the
output section of `MX17_README.md` for `p2` mode; the inherited MX17 modes
still write what they always did, plus the new generic per-volume scoring.

Each worker thread writes `<out>_t<N>.root` with three trees (or, without
ROOT, `<out>_t<N>_{events,clusters,volumes}.csv`). `hadd` them, or chain.

---

## 1. `EventTree` — one row per event

Unless `--skip-empty` is given, in which case only events with at least one
gas cluster are written. **The thrown count is printed at end of run and is
the correct normalization denominator** — the tree does not contain it.

### Core

| branch | unit | meaning |
|---|---|---|
| `eventID` | | |
| `edepDrift`, `edepAmp` | eV | energy deposited in `DriftGas` / `AmpGas` |
| `nPrimDrift`, `nPrimAmp` | | ion pairs from edep/W (see caveat §4) |
| `nClusDrift`, `nClusAmp` | | number of `ClusterTree` rows in each gas |
| `primInDrift`, `primInAmp` | bool | the primary track itself entered the gas |
| `primaryPDG`, `primaryCharge`, `primaryKE_MeV` | | the thrown particle |

### Classification — the branch the headline plot cuts on

| branch | meaning |
|---|---|
| `interactionClass` | integer code, see table below |
| `convProcess` | origin process of the dominant ionizing chain (`phot`, `compt`, `primary`, …) |
| `convVolume` | volume that chain started in |
| `convZ` | z of that interaction [mm] |
| `convEnergy` | KE of the particle born there [MeV] |
| `convParentEnergy` | creation KE of the *neutral* that converted [MeV] |
| `convParentBirthVolume` | where that neutral was born |

`interactionClass` codes (`InteractionClass` in `include/EventData.hh`):

| code | name | meaning |
|---|---|---|
| 0 | none | nothing deposited in the active gas |
| 1 | primary-track | a charged primary crossed the gas — **the MIP signal** |
| 2 | gas-photoelectric | photoabsorption inside `DriftGas`/`AmpGas` |
| 3 | gas-compton | Compton scatter inside the active gas |
| 4 | gas-other | pair production / other, in gas |
| 5 | wall-photoelectric | photoabsorption in a solid, electron reached the gas |
| 6 | wall-compton | Compton in a solid, electron reached the gas |
| 7 | wall-other | |
| 8 | fluorescence-reabs | soft (<20 keV) photon born in a solid, re-absorbed in gas — the Cu-K 8.05 / Fe-K 6.4 / Cr-K 5.4 keV chain |
| 9 | other | |

The label is assigned from whichever *ancestor chain* deposited the most
energy in the active gas. It is a convenience: every raw field it was
derived from is written too, so any finer cut can be made in analysis.

### First interaction of the primary

`firstIntProcess`, `firstIntVolume`, `firstIntX/Y/Z` [mm], `firstIntEnergy`
[MeV] — the first non-transport, non-msc, non-Rayleigh process the primary
underwent. For photon runs this is the conversion-layer budget at source,
independent of whether anything reached the gas.

### Per-layer sums [eV]

`edepFrontGas`, `edepCathGas`, `edepBackGas`, `edepWindowGas` — **gas in
field-free regions: deposits here produce no signal at all.** Recording
them is what separates "converted somewhere useless" from "did not convert".

`edepCathMylar`, `edepCathAl`, `edepWindowMylar`, `edepMesh`, `edepPadCu`,
`edepFR4`, `edepCuB`, `edepFrame` — the solids.

## 2. `ClusterTree` — one row per ionizing step in the active gas

`eventID`, `trackID`, `parentID`, `x/y/z` [mm], `edep` [eV], `nPrimary`,
`ke` [MeV], `volume`, `particle`, `time` [ns].

Provenance (added 2026-08-05):

| branch | meaning |
|---|---|
| `creatorProcess` | the process that made **this** track (`eIoni` for a delta ray) |
| `originProcess` | the process that started the **chain** — this is the one you want |
| `originVolume` | layer the chain started in |
| `originX/Y/Z` | where [mm] |
| `originEnergy` | KE of the particle born there [MeV] |
| `ancestorID` | trackID of that particle — group clusters by this |
| `parentEnergy` | creation KE of the converting neutral [MeV] |
| `parentBirthVolume` | where that neutral was born |

The distinction matters: a delta ray from a photoelectron has
`creatorProcess == "eIoni"` and a `parentID` pointing at the photoelectron,
so one-level ancestry loses the link to the conversion. `originProcess` /
`originVolume` are propagated down the whole track tree by
`src/TrackingAction.cc`. The rule: a track whose parent is *charged*
inherits its parent's origin; a track whose parent is *neutral* starts a new
origin. So "origin" = the interaction at which a neutral ancestor became
charged, and for electron runs it stays `primary` all the way down.

`parentEnergy` + `parentBirthVolume` are what identify the fluorescence
chain: a gas photoabsorption whose parent photon was born at 8.01 keV inside
`PCB_Cu_F` is Cu-K re-absorption, not a direct conversion of the beam photon.

## 3. `VolumeTree` — one row per (event, volume) with non-zero edep

`eventID`, `volume`, `edep` [eV], `nSteps`. Every volume, named or not.
Nothing is dropped. Use it to build the full energy budget of an event and
to catch deposits in layers nobody thought to add a named branch for.

## 4. Caveats that affect how these numbers should be read

- **`nPrimary` is `floor(edep/W)` plus a uniform-random remainder**, with
  one W per gas. That gives the right mean but **Bernoulli, not Fano,
  fluctuations**, and W is not constant below ~1 keV. It also ignores
  **Penning transfer**, which lowers the effective W in neon mixtures — so
  comparing Ar and Ne on `nPrimary` at fixed edep is biased. Prefer `edep`
  as the primary observable until Stage B does the conversion properly.
  (Campaign plan §5 step 3a.)
- **`edepAmp` is not proportional to signal.** An electron born at depth z
  in the 150 µm amplification gap is amplified only over the remaining
  distance, so the gain spans ~4 decades within that one volume. Use
  cluster `z` relative to the mesh; do not sum `edepAmp` and treat it like
  `edepDrift`.
- **Rayleigh scattering deposits nothing** and is excluded from
  `firstIntProcess`. Do not count it as an interaction.
- With `--skip-empty`, divide all rates by the **thrown** count from the run
  summary, not by `GetEntries()`.

## 5. Making the plot this was all built for

Event-by-event energy deposition, separated by interaction type:

```cpp
TChain c("EventTree"); c.Add("g60_ar_t*.root");
c.Draw("edepDrift+edepAmp >> hWall(100,0,20000)", "interactionClass==5");  // wall photoelectric
c.Draw("edepDrift+edepAmp >> hGasPE(100,0,20000)", "interactionClass==2", "same");
c.Draw("edepDrift+edepAmp >> hGasC(100,0,20000)",  "interactionClass==3", "same");
c.Draw("edepDrift+edepAmp >> hFluo(100,0,20000)",  "interactionClass==8", "same");
// ...and the MIP reference from an electron run:
TChain e("EventTree"); e.Add("e119_ar_t*.root");
e.Draw("edepDrift+edepAmp >> hMip(100,0,20000)", "interactionClass==1");
```

Conversion-layer budget (which layer the fakes come from):

```cpp
c.Draw("convVolume", "interactionClass>0");
```

Fluorescence lines, as a check that deexcitation is on:

```cpp
TChain cl("ClusterTree"); cl.Add("g60_ar_t*.root");
cl.Draw("parentEnergy*1000 >> (200,0,20)",
        "originProcess==\"phot\" && parentEnergy>0 && parentEnergy<0.02");
// expect Cu-K 8.05, Fe-K 6.40, Cr-K 5.41 keV
```
