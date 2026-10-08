// export_gas_json.C -- Magboltz .gas tables -> one JSON file for Stage B.
//
// Stage B is pure Python and must not need Garfield, so the transport
// coefficients are sampled here, on a log grid in E, and written as plain
// arrays. stage_b.gas.MagboltzTable interpolates them at the drift field.
//
//   source scripts/setup_garfield.sh
//   root -l -b -q 'scripts/garfield/export_gas_json.C("<dir with .gas>", "out.json")'
//
// Units in the JSON: E [V/cm], v [um/ns], D_L and D_T [um/sqrt(cm)],
// alpha and eta [1/cm].
R__LOAD_LIBRARY(libGarfield)
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "Garfield/MediumMagboltz.hh"

void export_gas_json(const char* dir = ".", const char* out = "magboltz.json") {
  const std::vector<std::string> gases = {"ArIso", "ArIso9010", "NeIso",
                                          "NeIso8515", "ArCO2Iso9352", "ArCF4Iso"};
  // Inside the table range only (100 V/cm - 90 kV/cm): no extrapolation.
  std::vector<double> E;
  for (int i = 0; i < 80; ++i) E.push_back(100. * std::pow(900., i / 79.));

  std::ofstream js(out);
  js << "{\n";
  for (size_t k = 0; k < gases.size(); ++k) {
    Garfield::MediumMagboltz m;
    const std::string file = std::string(dir) + "/" + gases[k] + ".gas";
    if (!m.LoadGasFile(file)) { std::cerr << "cannot load " << file << "\n"; return; }
    std::vector<double> v, dl, dt, al, et;
    for (double e : E) {
      double vx, vy, vz, l, t, a, h;
      m.ElectronVelocity(0, 0, -e, 0, 0, 0, vx, vy, vz);
      m.ElectronDiffusion(0, 0, -e, 0, 0, 0, l, t);
      m.ElectronTownsend(0, 0, -e, 0, 0, 0, a);
      m.ElectronAttachment(0, 0, -e, 0, 0, 0, h);
      v.push_back(std::fabs(vz) * 1e4);   // cm/ns -> um/ns
      dl.push_back(l * 1e4);              // sqrt(cm) -> um/sqrt(cm)
      dt.push_back(t * 1e4);
      al.push_back(a);
      et.push_back(h);
    }
    auto arr = [&](const char* name, const std::vector<double>& x, bool last = false) {
      js << "    \"" << name << "\": [";
      for (size_t i = 0; i < x.size(); ++i) js << (i ? ", " : "") << x[i];
      js << "]" << (last ? "\n" : ",\n");
    };
    js << "  \"" << gases[k] << "\": {\n";
    js << "    \"source\": \"" << gases[k] << ".gas (Magboltz, dry, 293.15 K, 760 Torr)\",\n";
    arr("E_V_per_cm", E);
    arr("v_um_per_ns", v);
    arr("DL_um_sqrtcm", dl);
    arr("DT_um_sqrtcm", dt);
    arr("alpha_per_cm", al);
    arr("eta_per_cm", et, true);
    js << "  }" << (k + 1 < gases.size() ? ",\n" : "\n");
  }
  js << "}\n";
  std::cout << "wrote " << out << std::endl;
}
