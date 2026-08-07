# Stage B — ionization clusters → per-pad charge and time

**Status (2026-08-07): skeleton complete and validated end-to-end.** It runs
on real Stage A output and produces per-pad charge and arrival time. Four
inputs are still placeholders, all listed below and all recorded in every
output file. **Do not quote a number from this stage without reading §3.**

```bash
python3 -m stage_b.selftest                       # 16 checks, no Geant4 needed
python3 -m stage_b.run stageA.root -o padhits.root
```

## 1. What it does

```
ClusterTree (Stage A)
  └─ depth from the drift frame        RunMeta.meshZ_mm − z_world
      ├─ attachment                    exp(−depth/λ)
      ├─ drift time                    depth / v_d
      ├─ longitudinal diffusion        N(0, σ_L √depth / v_d)   on time
      ├─ transverse diffusion          N(0, σ_T √depth)         on x, y
      ├─ mesh transparency             Bernoulli(ε)
      ├─ Polya gain                    Γ(θ+1, G/(θ+1))
      └─ pad lookup                    (x,y) → pad, polar
  └─ sum charge / earliest time per (event, pad)  →  PadHitTree
```

| module | role |
|---|---|
| `frame.py` | the drift frame, read from Stage A's `RunMeta` |
| `gas.py` | transport parameters; `TransportTable` is the Magboltz drop-in point |
| `padplane.py` | point → pad, built from `design/mapping/` |
| `clusters.py` | Stage A `ClusterTree` reader |
| `digitize.py` | the chain above |
| `run.py` | CLI; writes `PadHitTree` + `StageBMeta` |
| `selftest.py` | closed-form checks with synthetic clusters |

## 2. First result

400 muon events, Ar/iso, 3 mm, gain 10⁴, `--beam-spread 11.43`:

- **1.070 pads per event** (373 single / 26 double / 1 triple)
- median pad charge **37.6 fC**
- drift-time span **0–60 ns**, matching the 60 ns placeholder for Ar/iso

The pad multiplicity is the interesting one. **~1.1 pads/hit was previously
derived from a separate geometric analysis** (`vmm/nl_map.py`,
`SIM_CAMPAIGN_PLAN.md` §2); the full Geant4 → transport → pad chain now
reproduces it independently. That figure is what the "neighbour logic off"
baseline rests on, so having two independent routes to it matters.

## 3. Placeholders — what is not yet real

In rough order of how much each could move an answer.

1. **Gas transport parameters (P0.10).** `v_d`, σ_L, σ_T and attachment are
   placeholders from the campaign docs, shared with `vmm/time_resolution.py`
   so the two stay consistent. Anything scaling as 1/v_d — most of the timing
   — inherits their error. Swap in a `TransportTable` implementation when
   Magboltz tables exist; nothing else changes.
   *And expect the real dry table to still sit above the real detector:*
   MX17's bench finds **water contamination dominates v_d** (36.6 µm/ns
   measured vs a much higher dry prediction, 1–2 % H₂O inferred at SPS).
   `--v-scale` exists so a run states which assumption it was made under
   rather than hiding a preference in a constant.

2. **Induced current is not modelled (P0.17).** Each electron contributes its
   charge as a delta at its arrival time. The real signal is a fast electron
   spike plus a ~150 ns ion tail, and that shape sets the leading edge. This
   is the main model risk in the σ_t prediction, and it is cheap to fix for
   our non-resistive stack (static Ramo, Riegler closed form, already in
   Garfield++ `ComponentParallelPlate`).

3. **No pillars (P0.16).** They are not in the Geant4 geometry either, so
   electrons that should land on a pillar and be lost are amplified normally.
   Biases efficiency high by roughly the pillar area fraction, **~4.8 %**.

4. **Mesh transparency is a constant (P0.13).** Not a field map. Default 1.0,
   i.e. off; set `--mesh-transparency` explicitly if you want it.

Also absent by design, because Stage C owns them: shaping, threshold, ENC,
PDO/TDO. Use `vmm/` for those.

## 4. Traps

**World z is not drift depth.** The drift gas sits at positive world z and
electrons drift toward **+z**, so treating a cluster's z as "distance
drifted" is wrong by an offset *and* a sign — and stays finite and plausible
while being wrong. Stage A therefore records the frame
(`RunMeta.meshZ_mm` and friends) and `frame.py` converts. Never hardcode a z:
the drift gap is a scan axis. MX17 hit this and had half their gap at
negative depth.

**`nPrimary`, not `edep/W`.** Stage A already applied W and carried the
sub-W remainder probabilistically. Re-deriving would double-count.

**AmpGas clusters are not drift clusters.** They are already past the mesh,
so they skip the drift and are amplified over only the remaining gap:
`G^((gap−d)/gap)`. The gain spans ~4 decades across that one 150 µm volume.

**A pencil-beam Stage A file cannot give you pad multiplicity.** With
`beamSpread_mm = 0` every event lands on the same point, so the answer
describes that point, not the pad cell. `run.py` prints a NOTE when it sees
one. Stage A's default aim point is a pad centre in both polar coordinates —
which took two attempts, because fixing the radial coordinate put it in an
azimuthal gap.

**Mapping revision.** Pad *positions* are revision-independent (identical to
2.5 µm), so the geometry here is safe. Anything channel-level is not: the two
revisions agree on 11/1280 channel assignments. The revision used is recorded
in `StageBMeta`.

**Pad geometry source.** `padplane.py` reads `design/mapping/`, whose ring
radii sit 65–81 µm outside the gerber-derived pad mid-radii in
`include/P2PadMap.hh` (0.6 % of a pad). Irrelevant except within ~80 µm of an
edge; worth unifying on the gerber table when someone touches this next.
