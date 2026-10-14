// mesh_gain.cc -- avalanche gain of a COMSOL mesh cell (geometry study,
// 2026-10-14): the Townsend coefficient integrated along drift lines through
// the real mesh + amplification-gap field (DriftLineRKF::GetGain), with
// Penning transfer applied to the Magboltz table. MediumGas::
// EnablePenningTransfer(r, 0) rescales alpha with the excitation rates stored
// in the .gas file, so one table serves every transfer rate r.
//
// The gain is integrated here, not with DriftLineRKF::GetGain: the Magboltz
// table has no points between 3 and 15 kV/cm (drift / amplification grids,
// make_gas_table.C) and its interpolated alpha there is garbage (1e33 /cm at
// 10 kV/cm). alpha is taken from the table at E >= the first amplification
// grid point (15 kV/cm, alpha ~ 40 /cm) and 0 below; the few um of 3-15 kV/cm
// just under the mesh change ln G by < 0.01.
//
// Electrons start on an n x n grid over the cell, 20 um above the mesh; only
// lines that reach the pads (z < 5 um) count. Per r it prints
//   GAIN gas=... cell=... r=... n=... reached=... mean_gain=... mean_lnG=...
//        sd_lnG=... Eamp=...
// (Eamp: field at the hole centre, halfway down the gap, V/cm). Mean gain is
// the average over start points (what a track sees); no avalanche
// fluctuations, no space charge, no ion feedback.
//
// The measured gain slope (SPS DREAM mesh scans) fixes ln G; the gap and r are
// degenerate in the slope alone and are separated with the absolute gain.
//
//   mesh_gain <cell dir> <gas>.gas <out> <r1,r2,...> [n]
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Garfield/ComponentComsol.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"

using namespace Garfield;

namespace {
std::string Basename(const std::string& p) {
  std::string b = p.substr(p.find_last_of('/') + 1);
  return b.substr(0, b.rfind('.'));
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cerr << "usage: " << argv[0] << " <cell dir> <gas>.gas <out> <r1,r2,...> [n]\n";
    return 1;
  }
  const std::string dir = argv[1], gasfile = argv[2], out = argv[3];
  const int n = argc > 5 ? std::atoi(argv[5]) : 20;
  std::vector<double> rs;
  {
    std::stringstream ss(argv[4]);
    for (std::string t; std::getline(ss, t, ',');) rs.push_back(std::atof(t.c_str()));
  }

  // Materials as in mesh_transparency.cc (gas eps 1 in domains 1-3, eps 4 in 4-11).
  const std::string mpl = out + "_dielectrics.dat";
  {
    std::ofstream f(mpl);
    f << "2\n1 4\n11\n";
    for (int d = 1; d <= 11; ++d) f << d << " " << (d <= 3 ? 0 : 1) << "\n";
  }
  ComponentComsol fm;
  if (!fm.Initialise(dir + "/mesh.mphtxt", mpl, dir + "/potential.txt", "mm")) return 1;
  fm.EnablePeriodicityX();
  fm.EnablePeriodicityY();
  fm.EnableInterfaceCrossings(true);  // after the periodicities (patched build)

  // the bounding box is infinite with periodicity: use the unit cell
  double xmin, ymin, zmin, xmax, ymax, zmax;
  fm.GetElementaryCell(xmin, ymin, zmin, xmax, ymax, zmax);
  const double half = 0.5 * (xmax - xmin), xc = 0.5 * half;

  MediumMagboltz probe;  // only to find the mesh: any driftable medium
  probe.SetComposition("ar", 100.);
  fm.SetGas(&probe);
  double zlo = -1., zhi = -1.;
  for (double z = 0.; z <= 0.05; z += 0.00002) {
    Medium* med = fm.GetMedium(xc, xc, z);
    if (med == nullptr || !med->IsDriftable()) {
      if (zlo < 0.) zlo = z;
      zhi = z;
    }
  }
  if (zlo < 0.) {
    std::cerr << "mesh_gain: no mesh found\n";
    return 1;
  }
  const double z0 = zhi + 0.0020;
  double ex, ey, ez, v;
  Medium* med = nullptr;
  int st = 0;
  fm.ElectricField(0., 0., 0.5 * zlo, ex, ey, ez, v, med, st);
  const double eamp = std::sqrt(ex * ex + ey * ey + ez * ez);
  std::printf("cell half %.2f um, mesh %.1f - %.1f um, start %.1f um, "
              "Eamp(hole, mid-gap) %.0f V/cm\n", half * 1e4, zlo * 1e4,
              zhi * 1e4, z0 * 1e4, eamp);

  std::ofstream f(out + "_gain.txt");
  f << "# r x0 y0 [um] status z_end [um] gain\n";
  for (const double r : rs) {
    MediumMagboltz gas;
    if (!gas.LoadGasFile(gasfile)) return 1;
    if (r > 0. && !gas.EnablePenningTransfer(r, 0.)) {
      std::cerr << "mesh_gain: Penning transfer failed (no excitation rates?)\n";
      return 1;
    }
    fm.SetGas(&gas);
    Sensor sensor;
    sensor.AddComponent(&fm);
    sensor.SetArea(-10., -10., zmin, 10., 10., z0 + 0.01);
    DriftLineRKF drift(&sensor);
    drift.SetMaximumStepSize(2.e-5);
    drift.EnableAvalanche(false);
    std::vector<double> efs, bfs, afs;
    gas.GetFieldGrid(efs, bfs, afs);
    double eLow = 1.e9;
    for (double e : efs)
      if (e > 1.e4 && e < eLow) eLow = e;
    auto alpha = [&](double x, double y, double z) {
      double ex, ey, ez, v, a = 0.;
      Medium* md = nullptr;
      int s = 0;
      fm.ElectricField(x, y, z, ex, ey, ez, v, md, s);
      if (s != 0 || std::sqrt(ex * ex + ey * ey + ez * ez) < eLow) return 0.;
      gas.ElectronTownsend(ex, ey, ez, 0., 0., 0., a);
      return a;
    };
    int reached = 0;
    double sumG = 0., sumL = 0., sumL2 = 0.;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        const double x = -half + (i + 0.5) * 2. * half / n;
        const double y = -half + (j + 0.5) * 2. * half / n;
        drift.DriftElectron(x, y, z0, 0.);
        double x1, y1, z1, t1;
        int s1;
        drift.GetEndPoint(x1, y1, z1, t1, s1);
        // trapezoid on the drift-line points (steps <= 0.2 um)
        double lnG = 0., xa, ya, za, ta;
        drift.GetDriftLinePoint(0, xa, ya, za, ta);
        double aa = alpha(xa, ya, za);
        for (size_t k = 1; k < drift.GetNumberOfDriftLinePoints(); ++k) {
          double xb, yb, zb, tb;
          drift.GetDriftLinePoint(k, xb, yb, zb, tb);
          const double ab = alpha(xb, yb, zb);
          const double dl = std::sqrt((xb - xa) * (xb - xa) + (yb - ya) * (yb - ya) +
                                      (zb - za) * (zb - za));
          lnG += 0.5 * (aa + ab) * dl;
          xa = xb; ya = yb; za = zb; aa = ab;
        }
        const double g = std::exp(lnG);
        f << r << " " << x * 1e4 << " " << y * 1e4 << " " << s1 << " "
          << z1 * 1e4 << " " << g << "\n";
        if (z1 > 0.0005 || g <= 0.) continue;
        ++reached;
        sumG += g;
        sumL += lnG;
        sumL2 += lnG * lnG;
      }
    }
    const double m = reached ? sumL / reached : 0.;
    std::printf("GAIN gas=%s cell=%s r=%.2f n=%d reached=%d mean_gain=%.4g "
                "mean_lnG=%.4f sd_lnG=%.4f Eamp=%.0f\n", Basename(gasfile).c_str(),
                dir.c_str(), r, n * n, reached, reached ? sumG / reached : 0., m,
                reached ? std::sqrt(std::max(0., sumL2 / reached - m * m)) : 0., eamp);
  }
  return 0;
}
