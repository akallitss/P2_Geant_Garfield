// make_gas_table.C -- Magboltz gas table for one P2 mixture (task P0.10).
//
// Writes <out>.gas: drift velocity, longitudinal/transverse diffusion,
// Townsend and attachment coefficients on an E-field grid that covers both
// the drift gap (~100 V/cm - 3 kV/cm) and the 150 um amplification gap
// (~20 - 80 kV/cm). Stage B reads the drift part; the gain study reads the
// Townsend part. Penning transfer is NOT applied to the table here -- it only
// enters microscopic avalanches and is set there, per mixture, so it can be
// bracketed against the measured gain curves.
//
//   source scripts/setup_garfield.sh
//   root -l -b -q 'scripts/garfield/make_gas_table.C("ArIso9010", "out/ArIso9010", 10)'
//
// Gas names are the mm_sim names (src/GasMixtures.cc); the compositions below
// MUST match that table.
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libGarfield)
#include "Garfield/MediumMagboltz.hh"

struct Mix { std::vector<std::pair<std::string, double>> parts; };

void make_gas_table(const char* gas = "ArIso", const char* out = "ArIso",
                    int ncoll = 10, int nDrift = 24, int nAmp = 16) {
  std::map<std::string, Mix> mixtures = {
    {"ArIso",        {{{"ar", 95.}, {"ic4h10", 5.}}}},
    {"ArIso9010",    {{{"ar", 90.}, {"ic4h10", 10.}}}},
    {"NeIso",        {{{"ne", 95.}, {"ic4h10", 5.}}}},
    {"NeIso8515",    {{{"ne", 85.}, {"ic4h10", 15.}}}},
    {"ArCO2Iso9352", {{{"ar", 93.}, {"co2", 5.}, {"ic4h10", 2.}}}},
    {"ArCF4Iso",     {{{"ar", 88.}, {"cf4", 10.}, {"ic4h10", 2.}}}},
  };
  // Contaminated SPS gas (2026-10-14): Ar/CO2/iso 93/5/2 + O2 and/or H2O,
  // to test whether impurities explain the SPS drift scan. The base gas is
  // scaled by (1 - impurity) so the fractions add up to 100 %.
  // No O2/H2O monitor at SPS; 4 l/h through plastic tubing.
  auto contaminated = [](double ppmO2, double ppmH2O) {
    const double f = 1. - 1.e-6 * (ppmO2 + ppmH2O);
    Mix x{{{"ar", 93. * f}, {"co2", 5. * f}, {"ic4h10", 2. * f}}};
    if (ppmO2 > 0.) x.parts.push_back({"o2", 1.e-4 * ppmO2});
    if (ppmH2O > 0.) x.parts.push_back({"h2o", 1.e-4 * ppmH2O});
    return x;
  };
  mixtures["ArCO2Iso9352_O2_100"] = contaminated(100., 0.);
  mixtures["ArCO2Iso9352_O2_300"] = contaminated(300., 0.);
  mixtures["ArCO2Iso9352_O2_1000"] = contaminated(1000., 0.);
  mixtures["ArCO2Iso9352_H2O_500"] = contaminated(0., 500.);
  mixtures["ArCO2Iso9352_H2O_2000"] = contaminated(0., 2000.);
  mixtures["ArCO2Iso9352_O2_300_H2O_1000"] = contaminated(300., 1000.);
  auto it = mixtures.find(gas);
  if (it == mixtures.end()) {
    std::cerr << "make_gas_table: unknown gas " << gas << "\n";
    return;
  }
  const auto& p = it->second.parts;

  Garfield::MediumMagboltz m;
  auto g = [&](size_t i) { return i < p.size() ? p[i].first : std::string(""); };
  auto f = [&](size_t i) { return i < p.size() ? p[i].second : 0.; };
  m.SetComposition(g(0), f(0), g(1), f(1), g(2), f(2), g(3), f(3), g(4), f(4));
  m.SetTemperature(293.15);   // K
  m.SetPressure(760.);        // Torr; SPS/lab chambers run at ~1 atm

  // Field grid [V/cm]: log-spaced drift region + log-spaced amplification
  // region, joined into one sorted list.
  std::vector<double> efields;
  for (int i = 0; i < nDrift; ++i)
    efields.push_back(100. * std::pow(3000. / 100., double(i) / (nDrift - 1)));
  for (int i = 0; i < nAmp; ++i)
    efields.push_back(15000. * std::pow(90000. / 15000., double(i) / (nAmp - 1)));
  const std::vector<double> bfields = {0.};
  const std::vector<double> angles = {0.};
  m.SetFieldGrid(efields, bfields, angles);

  std::cout << "make_gas_table: " << gas << "  ncoll " << ncoll
            << "  " << efields.size() << " E points" << std::endl;
  m.GenerateGasTable(ncoll);
  m.WriteGasFile(std::string(out) + ".gas");
  std::cout << "wrote " << out << ".gas" << std::endl;
}
