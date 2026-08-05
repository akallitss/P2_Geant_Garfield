# Reference source files — provenance

Fetched 2026-08-05 from the public ATLAS Athena repository (Apache-2.0
license, © CERN for the benefit of the ATLAS collaboration), branch `main`,
path `MuonSpectrometer/MuonDigitization/MM_Digitization/`:

| file | role |
|---|---|
| `VMM_Shaper.{h,cxx}` | closed-form VMM shaper response, peak finding, threshold crossing, ballistic-deficit scaling (constants by G. Iakovidis, thesis CERN-THESIS-2014-148 §7.1.3) |
| `MM_ElectronicsResponseSimulation.{h,cxx}` | per-strip electronics loop: neighbor logic, time-window acceptance, dead-time-extended search window |
| `MM_DigitizationTool.h` | Gaudi property defaults (thresholds, dead times, cross-talk, NL flag) — values also summarized in `../README.md` |

These are kept verbatim as the ground truth for the Python port in
`../vmm_shaper.py`. Do not edit them.
