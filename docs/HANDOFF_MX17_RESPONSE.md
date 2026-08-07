# Handoff / feedback from the MX17 response-simulation effort

**Living document** — edited as the MX17 work progresses. Last update: **2026-08-07** (P2 worked through §2; see §2.5 for adoption status and §3 for asks going back).
Context: the MX17 detector (resistive-strip bulk Micromegas, DREAM readout) is building a full response-simulation chain, deliberately modeled on this repo's staged A/B/C architecture. Plan lives at `~/CLionProjects/MX17_Geant/design/RESPONSE_SIM_PLAN.md`. This file collects what flows *back* to the P2 campaign: suggestions, shared components, and review notes.

## 1. What MX17 adopted from P2 (so interfaces stay compatible)

- The staged, file-joined Stage A/B/C pipeline and its bookkeeping discipline (manifests, git+geometry hashes, no un-manifested runs).
- The Stage A `ClusterTree` schema (time + provenance + detector-frame coordinates) — MX17 is porting it verbatim. If the schema evolves here, tell MX17; if MX17 finds gaps, they'll be noted in this file.
- Tiered fidelity: microscopic Garfield++ once per (gas, HV) → frozen Polya/transparency parameters → fast sampler.
- The test-beam tuning firewall (`TESTBEAM_PLAN.md` §4–6): minimal global parameter set, blind predictions first, named escalation paths instead of per-plot tuning.

## 2. Suggestions for the P2 campaign (current)

### 2.1 Induction upgrade is cheaper than Phase 4 assumes
`TIME_RESOLUTION_NOTES.md` flags the delta-charge-into-shaper model as the main σ_t risk and defers "realistic induced current (electron peak + ~150 ns ion tail)" to Phase 4. MX17 is building a semi-spectral quasi-static weighting-potential solver for its resistive stack (see plan §3). For P2's *non-resistive* stack this machinery degenerates to a **static** N-layer weighting potential — no time dependence, no FEM: Riegler's closed-form layered solutions (JINST 11 (2016) P11002), already implemented in Garfield++ `ComponentParallelPlate` (`AddPixel`/`AddStrip` on an N-layer stack). Electron spike + analytic uniform-field ion current through that weighting potential gives the realistic i(t) at roughly the cost of one afternoon, not a phase. Concretely for Stage B: replace each electron's delta charge with `i(t) = Q_e·δ_fast(t) + Q_ion·i_ion(t; μ_ion, gap)` weighted by the pad's Ramo potential — the ion tail then interacts correctly with the VMM shaper instead of via the hardcoded 150 ns ballistic-deficit factor.
- Related nit: `MM_ION_FLOW_TIME_NS = 150.0` (Athena constant in `vmm/vmm_shaper.py`) encodes ATLAS NSW's gap/gas ion time. P2's gap and gas differ; once real induction exists, drop the constant rather than retune it.

### 2.2 Shared Stage B proposal
P2's Stage B is not yet written; MX17's is being written now (Python, numpy/uproot, file contracts in the MX17 plan §2, §7). The stage decomposition is identical up to two plug-ins: (a) induction kernel (static Ramo for P2, dynamic resistive templates for MX17), (b) electronics (VMM vs DREAM). Proposal: when P2 starts Stage B, lift `response/digitizer/` from MX17_Geant rather than writing fresh — drift/diffusion/attachment sampling, Polya gain, mesh transparency, packet bookkeeping, and the manifest tooling are detector-agnostic. This file will note when that package is stable enough to lift (not yet — created 2026-08-06).

### 2.3 Mesh unit-cell field map (S2) is directly reusable
MX17 is solving a woven-mesh unit cell (Garfield++ neBEM, desktop-scale job) for transparency ε(E_d/E_a), funneling, and ion endpoints — replacing the "constant ε_mesh in baseline, field map in Phase 4" placeholder. P2's mesh differs only in weave parameters; the same script with different constants covers the P2 Phase-4 item (`SIM_CAMPAIGN_PLAN.md:789`) and the measured-transparency escalation path in `TESTBEAM_PLAN.md` §6. Will link the script path here when it exists.

### 2.4 Smaller review notes on the current P2 code/docs
- `vmm/time_resolution.py` gas table `v_d`/`σ_L` placeholders: MX17's bench experience says **water contamination dominates** drift velocity in real gas systems (measured 36.6 µm/ns vs far higher dry-Magboltz prediction in Ar/iso; 1–2% H2O inferred at SPS). When the Magboltz tables get generated (P0.10), generate wet variants (0.5/1/2% H2O) alongside dry — cheap, and the bench data will land between them.
- ENC ask (`NEEDED_INPUTS.md` item 5): MX17's DREAM pedestal-run approach (derive per-channel σ AND the common-mode covariance from dedicated pedestal acquisitions, inject correlated noise, let the analysis CNS remove it as in data) transfers to VMM benches; worth requesting pedestal-style runs in the P2 test-beam shopping list (`TESTBEAM_PLAN.md` §2 item 1 already asks for pedestal σ — add the covariance/common-mode request).
- **Pillar material is Dynamask** (photoimageable dry-film solder mask) — per Dylan (2026-08-06), for the MX17 modules **and for the P2 wedge pillars too**. Relevant here for: the Geant4 pillar material definition (epoxy-acrylate dry film, ρ ≈ 1.2–1.4 g/cm³ — not kapton/FR4), the material-budget accounting, and pillar dead-spot modeling in Stage B (P2 gerber pillar layout × Dynamask dielectric). If a mesh-unit-cell field solve ever includes a pillar, use Dynamask ε_r ≈ 3.9.
- The edep→ADC waterfall audit (`SIM_CAMPAIGN_PLAN.md` §"proportionality") is excellent; MX17 adds one row for resistive detectors (RC spreading × per-channel threshold interaction). If P2 ever tests a resistive prototype, take the row back.

## 2.5 Adoption status of §2 (P2 side, updated 2026-08-07)

| suggestion | status in the P2 repo |
|---|---|
| §2.1 induction upgrade is cheap | ✅ **adopted.** Promoted out of Phase 4 to **P0.17**; `TIME_RESOLUTION_NOTES.md` §4.1 and §6 rewritten with the static-Ramo/Riegler argument, and `SIM_CAMPAIGN_PLAN.md` §8 now records *why* the original cost estimate was wrong for a non-resistive stack. The `MM_ION_FLOW_TIME_NS` nit is captured in P0.17 as "drop, don't retune". |
| §2.2 shared Stage B | ✅ **adopted as a standing instruction.** P0.9 now says lift `response/digitizer/` rather than write fresh, and to **re-check this file before starting** — the package was explicitly not stable as of 2026-08-06. |
| §2.3 mesh unit-cell field map | ✅ **adopted.** Recorded against the Phase-4 mesh-transparency bullet and cross-linked to P0.13 and `TESTBEAM_PLAN.md` §6, with "ask before writing our own". Waiting on MX17 for the script path. |
| §2.4 wet Magboltz variants | ✅ **adopted.** P0.10 now generates 0.5/1/2 % H₂O alongside dry; `TIME_RESOLUTION_NOTES.md` §4.2 and §6 say to quote σ_t as a dry→2 % bracket rather than a line. |
| §2.4 pedestal common-mode covariance | ✅ **adopted.** `TESTBEAM_PLAN.md` §2 item 1 now asks for dedicated pedestal runs **and** the channel-to-channel covariance, with the reason spelled out (independent noise averages down across a cluster, common-mode does not, so σ-only is optimistic exactly at low threshold). |
| §2.4 Dynamask pillar material | ✅ adopted earlier (2026-08-06): `HANDOFF.md` §3, `NEEDED_INPUTS.md` status table, P0.16. |
| §2.4 waterfall row for resistive detectors | n/a unless P2 tests a resistive prototype; noted, not filed. |

The three readout-accuracy upgrades from the separate MX17_Geant model work
(`NEEDED_INPUTS.md` "Readout-copper zoning…") were documented but never
ticketed; they are now **P0.18** (Cu zoning — a live error in P2, and the one
that matters most since copper is ~90 % of the argon photon fake budget),
**P0.19** (real pad artwork) and **P0.20** (resistive layer structure).

## 3. Asks from P2 to MX17

1. **§2.2 — say when `response/digitizer/` is stable enough to lift**, and
   note it in this file. P0.9 is now written to depend on that signal, so
   the failure mode if it never arrives is that we duplicate a week of
   detector-agnostic work. A rough date is more useful than silence.
2. **§2.3 — the mesh unit-cell neBEM script path**, plus which weave
   parameters are constants vs inputs. P2's mesh is 48/19 µm plain weave
   (confirmed 2026-08-05); if the script takes wire/pitch as arguments we
   get P0.13 and the Phase-4 transparency item nearly free.
3. **Stage A `ClusterTree` schema change, per §1's request to be told.**
   The schema MX17 is porting has since grown full provenance (creator
   process, origin process/volume/position/energy, ancestor trackID,
   converting parent's birth energy and volume). Details in
   `OUTPUT_FORMAT.md` §2.
4. **New this session: a `RunMeta` tree** (`OUTPUT_FORMAT.md` §0), one row
   per worker file, carrying thrown count, git hash + dirty flag, a
   geometry hash + readable digest, resolved gas composition and the beam
   configuration. This is the piece that makes "no un-manifested runs"
   enforceable rather than aspirational, and it is detector-agnostic —
   worth MX17 copying, and worth keeping compatible if it does.
5. **Which pure-gas W values MX17 uses, and their source.** Ours are
   inherited and unsourced (P0.2c); isobutane in particular is carried at
   26.0 eV against a commonly quoted 23.4 eV. If MX17 has already done this
   sourcing, we take it rather than repeat it.
6. **Water contamination (§2.4): how was the 1–2 % H₂O at SPS inferred**,
   and from what measurement? If it is a drift-velocity fit rather than a
   direct hygrometer reading, that changes how we bracket it — and whether
   the P2 test-beam shopping list should ask for gas-quality monitoring.

## 4. Change log
- 2026-08-06: created; initial suggestions §2.1–2.4 based on survey of this repo at commit `70bc5a8` + working tree of 2026-08-06.
- 2026-08-06 (later): added pillar-material note in §2.4 — Dynamask for both MX17 and P2 pillars (from Dylan).
- 2026-08-07 (P2 side): worked through §2. Added §2.5 adoption table; filed §2.1 as P0.17 and the three readout upgrades as P0.18–P0.20; folded §2.2/§2.3/§2.4 into P0.9, P0.10, Phase 4, `TIME_RESOLUTION_NOTES.md` and `TESTBEAM_PLAN.md`; filled §3 with six asks back.
