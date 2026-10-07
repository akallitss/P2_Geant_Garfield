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
  const std::map<std::string, Mix> mixtures = {
    {"ArIso",        {{{"ar", 95.}, {"ic4h10", 5.}}}},
    {"ArIso9010",    {{{"ar", 90.}, {"ic4h10", 10.}}}},
    {"NeIso",        {{{"ne", 95.}, {"ic4h10", 5.}}}},
    {"NeIso8515",    {{{"ne", 85.}, {"ic4h10", 15.}}}},
    {"ArCO2Iso9352", {{{"ar", 93.}, {"co2", 5.}, {"ic4h10", 2.}}}},
    {"ArCF4Iso",     {{{"ar", 88.}, {"cf4", 10.}, {"ic4h10", 2.}}}},
  };
  auto it = mixtures.find(gas);
  if (it == mixtures.end()) {
    std::cerr << "make_gas_table: unknown gas " << gas << "\n";
    return;
  }
  const auto& p = it->second.parts;

  Garfield::MediumMagboltz m;
  if (p.size() == 2)
    m.SetComposition(p[0].first, p[0].second, p[1].first, p[1].second);
  else
    m.SetComposition(p[0].first, p[0].second, p[1].first, p[1].second,
                     p[2].first, p[2].second);
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
