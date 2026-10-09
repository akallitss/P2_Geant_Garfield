// SolveP2MeshCell.java -- re-solve the woven-mesh unit cell for P2 and write
// the Garfield++ ComponentComsol inputs.
//
// Base model: MeshUnitCellIon_PressureStudy.mph (colleague's mesh cell,
// COMSOL 6.3; used with permission, 2026-10-09). Only parameters change:
//   pitch   2*68 um  (measured P2 mesh: 18 um wire, 50 um opening; the
//                     weave repeats every two wires, hence 2*)
//   rWire   9 um     (unchanged)
//   meshHeight 150 um (amplification gap, unchanged)
//   P2 is polarised on mesh and drift, pads grounded (SPS July 2026):
//   VAnode 0 V, VMesh -450 V, VCathode = -450 V - E_drift * height
//
// Two cases, same 625 V/cm drift field:
//   h4mm : height 4 mm, VCathode -700 V   (the real drift gap)
//   h2mm : height 2 mm, VCathode -575 V   (the original cell height)
//
// Per case it writes, into OUT/<case>/:
//   mesh.mphtxt               mesh1 export
//   potential.txt             data1 (V)
//   wpot_anode.txt            data2 (V2, anode weighting potential)
//   wpot_mesh.txt             data3 (V4, mesh weighting potential)
//   wpot_cathode.txt          data4 (V5, cathode weighting potential)
// and prints the domain -> permittivity table (MATERIAL lines) to the log.
//
//   comsol compile SolveP2MeshCell.java
//   comsol batch -np 4 -nosave -inputfile SolveP2MeshCell.class
import com.comsol.model.*;
import com.comsol.model.util.*;

public class SolveP2MeshCell {
  static final String IN =
      "/eos/user/a/akallits/p2_sim/comsol/MeshUnitCellIon_PressureStudy.mph";
  static final String OUT = "/eos/user/a/akallits/p2_sim/comsol/p2_sps";

  public static void main(String[] args) throws java.io.IOException { run(); }

  public static Model run() throws java.io.IOException {
    Model m;
    String[][] cases = {{"h4mm", "4[mm]", "-700[V]"}, {"h2mm", "2[mm]", "-575[V]"}};
    for (String[] c : cases) {
      m = ModelUtil.load("Model", IN);
      m.param().set("pitch", "2*68[um]");
      m.param().set("height", c[1]);
      m.param().set("VAnode", "0[V]");
      m.param().set("VMesh", "-450[V]");
      m.param().set("VCathode", c[2]);

      String dir = OUT + "/" + c[0];
      // the folder is made by run_mesh_cell.sh: COMSOL forbids Java file I/O
      System.out.println("== case " + c[0] + ": height " + c[1] + ", VCathode " + c[2]);
      long t0 = System.currentTimeMillis();
      m.study("std1").run();
      System.out.println("   std1 solved in " + (System.currentTimeMillis() - t0) / 1000 + " s");

      String[][] ex = {{"mesh1", "mesh.mphtxt"}, {"data1", "potential.txt"},
                       {"data2", "wpot_anode.txt"}, {"data3", "wpot_mesh.txt"},
                       {"data4", "wpot_cathode.txt"}};
      for (String[] e : ex) {
        m.result().export(e[0]).set("filename", dir + "/" + e[1]);
        m.result().export(e[0]).run();
        System.out.println("   wrote " + e[1]);
      }

      // Domain -> permittivity table for Garfield's ComponentComsol, printed
      // to the job log (COMSOL's security settings forbid writing files).
      for (String mt : m.component("comp1").material().tags()) {
        Material mat = m.component("comp1").material(mt);
        String eps;
        try { eps = String.join(" ", mat.propertyGroup("def").getStringArray("relpermittivity")); }
        catch (Exception ex2) { eps = "?"; }
        StringBuilder sb = new StringBuilder();
        for (int d : mat.selection().entities()) sb.append(d).append(' ');
        System.out.println("   MATERIAL " + mt + " | " + mat.label() + " | eps_r " + eps
                           + " | domains " + sb.toString().trim());
      }
      ModelUtil.remove("Model");
    }
    return null;   // nothing to save: run with -nosave
  }
}
