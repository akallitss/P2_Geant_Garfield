// mesh_transparency.C -- field map and electron transparency of the P2 woven
// mesh, from the COMSOL unit cell (scripts/comsol/, task P0.13).
//
// The cell is one weave period (2 x 68 um) in x and y, periodic, with the
// pads at z = 0, the mesh around z = 150 um and the drift electrode on top.
// The probe mode shows the wires at x = +-34 um (along y) and y = +-34 um
// (along x), solid between z = 135.5 and 165 um where two wires cross and
// 141-159.5 um on a single wire; the openings are centred on (0, 0).
// Free-flight steps that clip a wire are caught only with the field map's
// interface panels (ComponentFieldMap::EnableInterfaceCrossings, pinned build
// 927e5c21 or later; off by default). Runs before 2026-10-14 did not enable
// them, so only step end points inside a wire were caught.
//
// Built and run by run_transparency.sh (a compiled program: ROOT 6.38's
// cling fails to JIT this against Garfield):
//   mesh_transparency probe  <cell dir> "" <out>
//   mesh_transparency maps   <cell dir> <gas>.gas <out>
//   mesh_transparency transp <cell dir> <gas>.gas <out> [n] [seed] [flags]
// flags (comma-separated): nox  = no interface-crossing check (old runs),
//                          null = AvalancheMicroscopic null-collision steps
//
// probe : prints the medium along z at a few (x, y) to locate the wires.
// maps  : writes <out>_xz_y<y>.txt and <out>_xy.txt (potential, field, inside
//         gas or not) on regular grids, and <out>_lines.txt: field lines
//         (RKF drift lines without diffusion) started across the cell.
// transp: <n> electrons uniform in (x, y) at kStartAbove above the mesh,
//         tracked microscopically (with diffusion) until they cross below the
//         mesh or stop on a wire; writes one line per electron to
//         <out>_electrons.txt and prints the transparency. The same starting
//         points are also run with RKF (field lines only) for comparison.
//
// Lengths below are in cm (Garfield's unit); the COMSOL export is in mm.
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"

namespace {
const double kHalf = 0.0068;       // half cell [cm] (68 um)
// Above / below the woven wires [cm]: found from the map at a wire crossing
// (solid extent +- 5 um); 170 / 130 um for the colleague's meshHeight 150.
double kMeshTop = 0.0170;
double kMeshBottom = 0.0130;
const double kStartAbove = 0.015;  // start this far above kMeshTop [cm]

struct Mix { std::vector<std::pair<std::string, double>> parts; };
const std::map<std::string, Mix> kMixtures = {
  {"ArIso",        {{{"ar", 95.}, {"ic4h10", 5.}}}},
  {"ArIso9010",    {{{"ar", 90.}, {"ic4h10", 10.}}}},
  {"NeIso",        {{{"ne", 95.}, {"ic4h10", 5.}}}},
  {"NeIso8515",    {{{"ne", 85.}, {"ic4h10", 15.}}}},
  {"ArCO2Iso9352", {{{"ar", 93.}, {"co2", 5.}, {"ic4h10", 2.}}}},
  {"ArCF4Iso",     {{{"ar", 88.}, {"cf4", 10.}, {"ic4h10", 2.}}}},
  // Literature comparison only (Nikolopoulos et al., JINST 6 (2011) P06011);
  // no gas table: microscopic tracking only, RKF skipped.
  {"ArCO2_8515",   {{{"ar", 85.}, {"co2", 15.}}}},
};

std::string Basename(const std::string& p) {
  std::string b = p.substr(p.find_last_of('/') + 1);
  return b.substr(0, b.find('.'));
}
}  // namespace

void mesh_transparency(const char* mode = "probe", const char* cell = ".",
                       const char* gasfile = "", const char* out = "mt",
                       int n = 2000, int seed = 1, const char* flags = "") {
  using namespace Garfield;
  const std::string m(mode), dir(cell);

  // Gas: composition from the mm_sim name, transport table for RKF.
  MediumMagboltz gas;
  bool haveTable = false;
  gas.SetTemperature(293.15);
  gas.SetPressure(760.);
  if (std::string(gasfile) != "") {
    const std::string name = Basename(gasfile);
    auto it = kMixtures.find(name);
    if (it == kMixtures.end()) {
      std::cerr << "mesh_transparency: unknown gas " << name << "\n";
      return;
    }
    const auto& p = it->second.parts;
    if (p.size() == 2)
      gas.SetComposition(p[0].first, p[0].second, p[1].first, p[1].second);
    else
      gas.SetComposition(p[0].first, p[0].second, p[1].first, p[1].second,
                         p[2].first, p[2].second);
    haveTable = std::ifstream(gasfile).good() && gas.LoadGasFile(gasfile);
    if (!haveTable)
      std::cout << "mesh_transparency: no gas table " << gasfile
                << ", microscopic tracking only\n";
    gas.SetMaxElectronEnergy(200.);
    gas.Initialise(true);
  } else {
    gas.SetComposition("ar", 95., "ic4h10", 5.);
  }

  // Materials: written next to the output (gas eps 1 in domains 1-3,
  // mesh/insulator eps 4 in 4-11, as printed by SolveP2MeshCell).
  const std::string mpl = std::string(out) + "_dielectrics.dat";
  {
    std::ofstream f(mpl);
    f << "2\n1 4\n11\n";
    for (int d = 1; d <= 11; ++d) f << d << " " << (d <= 3 ? 0 : 1) << "\n";
  }
  ComponentComsol fm;
  if (!fm.Initialise(dir + "/mesh.mphtxt", mpl, dir + "/potential.txt", "mm"))
    return;
  fm.EnablePeriodicityX();
  fm.EnablePeriodicityY();
  const std::string fl(flags);
  const bool xings = fl.find("nox") == std::string::npos;
  if (xings) {
    // after the periodicities: the panels are built from the current setup
    fm.EnableInterfaceCrossings(true);
    std::cout << "interface panels: " << fm.GetNumberOfInterfacePanels() << "\n";
  }
  fm.SetGas(&gas);
  fm.PrintRange();
  fm.PrintMaterials();

  double xmin, ymin, zmin, xmax, ymax, zmax;
  fm.GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax);
  const double zTop = zmax;

  auto inGas = [&](double x, double y, double z) {
    Medium* med = fm.GetMedium(x, y, z);
    return med != nullptr && med->IsDriftable();
  };
  {
    double zlo = -1., zhi = -1.;
    for (double z = 0.; z <= 0.05; z += 0.00002) {
      if (!inGas(0.0034, 0.0034, z)) {
        if (zlo < 0.) zlo = z;
        zhi = z;
      }
    }
    if (zlo > 0.) {
      kMeshBottom = zlo - 0.0005;
      kMeshTop = zhi + 0.0005;
    }
    std::printf("mesh solid at the crossing: %.1f - %.1f um; counted as passed "
                "below %.1f um, start %.1f um\n", zlo * 1e4, zhi * 1e4,
                kMeshBottom * 1e4, (kMeshTop + kStartAbove) * 1e4);
  }

  if (m == "probe") {
    const double pts[][2] = {{0., 0.}, {0.0034, 0.0034}, {0.0034, 0.},
                             {0., 0.0034}, {-0.0034, 0.0034}, {0.005, 0.002}};
    for (auto& p : pts) {
      std::printf("x=%6.1f um y=%6.1f um : solid z [um] =", p[0] * 1e4,
                  p[1] * 1e4);
      bool prev = true;
      for (double z = 0.; z <= 0.03; z += 0.00005) {
        const bool g = inGas(p[0], p[1], z);
        if (g != prev) std::printf(" %s%.1f", g ? "..." : " ", z * 1e4);
        prev = g;
      }
      double ex, ey, ez, v;
      Medium* med = nullptr;
      int st = 0;
      fm.ElectricField(p[0], p[1], 0.1, ex, ey, ez, v, med, st);
      std::printf("   | E(z=1mm) = %.0f V/cm, V = %.1f\n", -ez, v);
    }
    for (double z : {0.005, 0.03, 0.1, 0.2, zTop - 0.01}) {
      double ex, ey, ez, v;
      Medium* med = nullptr;
      int st = 0;
      fm.ElectricField(0.0034, 0.0034, z, ex, ey, ez, v, med, st);
      std::printf("z=%7.1f um : Ez = %8.0f V/cm  V = %8.2f\n", z * 1e4, ez, v);
    }
    return;
  }

  Sensor sensor;
  sensor.AddComponent(&fm);

  if (m == "maps") {
    auto writeGrid = [&](const std::string& fn, bool xz, double fixed,
                         double a0, double a1, double b0, double b1, int na,
                         int nb) {
      std::ofstream f(fn);
      f << "# " << (xz ? "x z" : "x y") << " [um], V [V], Ex Ey Ez [V/cm], gas\n";
      for (int j = 0; j < nb; ++j) {
        const double b = b0 + (b1 - b0) * j / (nb - 1);
        for (int i = 0; i < na; ++i) {
          const double a = a0 + (a1 - a0) * i / (na - 1);
          const double x = a, y = xz ? fixed : b, z = xz ? b : fixed;
          double ex = 0, ey = 0, ez = 0, v = 0;
          Medium* med = nullptr;
          int st = 0;
          fm.ElectricField(x, y, z, ex, ey, ez, v, med, st);
          const int g = (st == 0 && med && med->IsDriftable()) ? 1 : 0;
          f << a * 1e4 << " " << b * 1e4 << " " << v << " " << ex << " " << ey
            << " " << ez << " " << g << "\n";
        }
      }
      std::cout << "wrote " << fn << "\n";
    };
    // Side views through the cell (y = 0 and y = 17, 34 um), pads to 400 um.
    for (double y : {0., 0.0017, 0.0034})
      writeGrid(std::string(out) + "_xz_y" + std::to_string(int(y * 1e4 + 0.5)) +
                    ".txt", true, y, -kHalf, kHalf, 0., 0.04, 273, 801);
    // Top views: through the wires, just above and just below them.
    for (double z : {0.0150, 0.0185, 0.0110})
      writeGrid(std::string(out) + "_xy_z" + std::to_string(int(z * 1e4 + 0.5)) +
                    ".txt", false, z, -kHalf, kHalf, -kHalf, kHalf, 273, 273);

    // Field lines: RKF drift lines without diffusion from 400 um, across x at
    // two y values; projected later onto the x-z plane.
    if (!haveTable) return;
    DriftLineRKF drift(&sensor);
    drift.SetMaximumStepSize(2.e-5);
    std::ofstream f(std::string(out) + "_lines.txt");
    f << "# line y0 x y z [um]\n";
    int line = 0;
    for (double y0 : {0., 0.0034}) {
      for (int i = 0; i < 41; ++i) {
        const double x0 = -kHalf + 2 * kHalf * (i + 0.5) / 41.;
        drift.DriftElectron(x0, y0, 0.04, 0.);
        const auto np = drift.GetNumberOfDriftLinePoints();
        for (size_t k = 0; k < np; ++k) {
          double x, y, z, t;
          drift.GetDriftLinePoint(k, x, y, z, t);
          f << line << " " << y0 * 1e4 << " " << x * 1e4 << " " << y * 1e4
            << " " << z * 1e4 << "\n";
        }
        ++line;
      }
    }
    std::cout << "wrote " << out << "_lines.txt (" << line << " lines)\n";
    return;
  }

  if (m == "transp") {
    Random::SetEngine(RandomEngineRoot(seed));
    const double z0 = kMeshTop + kStartAbove;
    // Stop once below the wires; the amplification gap is not tracked here.
    // The top is 1 mm above the start so back-diffusion cannot end a track.
    sensor.SetArea(-10., -10., kMeshBottom - 0.0005, 10., 10., z0 + 0.1);
    AvalancheMicroscopic aval(&sensor);
    if (fl.find("null") != std::string::npos) aval.EnableNullCollisionSteps(true);
    DriftLineRKF drift(&sensor);
    drift.SetMaximumStepSize(2.e-5);

    std::ofstream f(std::string(out) + "_electrons.txt");
    f << "# x0 y0 [um] | micro: status x1 y1 z1 [um] pass | rkf: status x1 y1 z1 [um] pass\n";
    int passM = 0, passR = 0, attached = 0, wire = 0, other = 0;
    for (int i = 0; i < n; ++i) {
      const double x = -kHalf + 2 * kHalf * RndmUniform();
      const double y = -kHalf + 2 * kHalf * RndmUniform();
      aval.DriftElectron(x, y, z0, 0., 0.1);
      double a0, b0, c0, t0, e0, x1, y1, z1, t1, e1;
      int st = 0;
      aval.GetElectronEndpoint(0, a0, b0, c0, t0, e0, x1, y1, z1, t1, e1, st);
      const int pm = z1 < kMeshBottom ? 1 : 0;
      passM += pm;
      if (!pm) {
        if (st == -7) ++attached;  // attachment, not the mesh
        else if (st == -5) ++wire; // stopped in a wire
        else ++other;              // e.g. left through the top
      }
      double xr = 0, yr = 0, zr = 0, tr = 0;
      int sr = 0;
      if (haveTable) {
        drift.DriftElectron(x, y, z0, 0.);
        drift.GetEndPoint(xr, yr, zr, tr, sr);
      }
      const int pr = zr < kMeshBottom ? 1 : 0;
      passR += pr;
      f << x * 1e4 << " " << y * 1e4 << " " << st << " " << x1 * 1e4 << " "
        << y1 * 1e4 << " " << z1 * 1e4 << " " << pm << " " << sr << " "
        << xr * 1e4 << " " << yr * 1e4 << " " << zr * 1e4 << " " << pr << "\n";
      if ((i + 1) % 200 == 0)
        std::printf("%d/%d  micro %.4f  rkf %.4f\n", i + 1, n,
                    double(passM) / (i + 1), double(passR) / (i + 1));
    }
    const double tm = double(passM) / n, tr = double(passR) / n;
    std::printf("TRANSPARENCY gas=%s cell=%s n=%d micro=%.4f +- %.4f rkf=%.4f "
                "attached=%.4f wire=%.4f other=%.4f flags=%s\n",
                Basename(gasfile).c_str(), dir.c_str(), n, tm,
                std::sqrt(tm * (1 - tm) / n), haveTable ? tr : -1.,
                double(attached) / n, double(wire) / n, double(other) / n,
                xings ? fl.c_str() : "nox");
    return;
  }
  std::cerr << "mesh_transparency: unknown mode " << m << "\n";
}

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cerr << "usage: " << argv[0]
              << " probe|maps|transp <cell dir> <gasfile|\"\"> <out> [n] [seed] [flags]\n";
    return 1;
  }
  mesh_transparency(argv[1], argv[2], argv[3], argv[4],
                    argc > 5 ? std::atoi(argv[5]) : 2000,
                    argc > 6 ? std::atoi(argv[6]) : 1, argc > 7 ? argv[7] : "");
  return 0;
}
