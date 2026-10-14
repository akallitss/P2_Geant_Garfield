# DREAM charge calibration and the measured absolute gain (2026-10-14)

## Charge range actually used
- DREAM gain is set per channel by registers 6/7 (2 bits per channel; DREAM
  user manual v3, tables 32-33): 00 = 50 fC, 01 = 100 fC, 10 = 200 fC, 11 = 600 fC.
- SPS July 2026 (`P2B_Beam.cfg`, e.g. runs/run_31/.../raw_daq_data) and the cosmic
  bench (`CosmicTb_P2.cfg`): `Dream * 6/7 = 0xAAAA` -> **200 fC on every channel**.
- Register 1 bits 2-3 ("100 fF") are the internal TEST capacitor, not the gain.

## ADC per fC
- Manual: 2 V p-p output, 12-bit ADC over +-1 V -> nominal 20 mV/fC at 100 fC.
- Measured (I. Mandjavidze, 220208_Dream2TpotTb.pptx, 100 fC range, 180 ns):
  18.0 mV/fC at 0 pF, 17.9 at 100 pF; noise 23 ADC = ENC 4000 e -> 37 ADC/fC.
- 200 fC range: **~18 ADC/fC, 1 ADC ~ 340 e** (+-10 %, delta input at the peak).

## Absolute effective gain (charge seen at the 180 ns peak)
| data | gas | mesh | signal | primaries | G_eff |
|---|---|---|---|---|---|
| Fe55, cosmic det2 (stage 18) | Ar/iso 95/5 | 415 V | photopeak 994 ADC (cluster sum) | 227 (5.9 keV / 26 eV) | ~1.5e3 |
| MIP, SPS P2_OUT (beam_nominal_meshscan_1) | Ar/CO2/iso 93/5/2 | 450 V | Landau MPV 330 ADC | 17-29 (Heed 4 mm MPV-median) x T 0.975 | 4.0-6.9e3 |

The ion signal is longer than the 180 ns peaking time: for ions crossing the
150 um gap in ~250-500 ns a CR-RC^4 peak collects 0.46-0.76 of the charge
(box-current estimate), so the true gain is ~1.3-2.2 x G_eff:
Fe55 415 V: G ~ 2-3e3 (ln G 7.6-8.1); SPS 450 V: G ~ 0.5-1.5e4 (ln G 8.6-9.6).

Consistency with the slope: ln G x dln(alpha)/dln(E) = V / V_e
(Fe55 415/28 = 14.8, SPS 450/25 = 18) -> dln(alpha)/dln(E) ~ 1.9-2.0, as
Magboltz gives at 27-30 kV/cm (2.1). The simulation (mesh_gain.cc) must
reach these ln G at the measured voltages with one Penning r and one gap.
