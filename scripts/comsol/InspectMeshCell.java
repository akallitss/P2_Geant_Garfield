// InspectMeshCell.java -- print what the mesh unit-cell model contains:
// parameters, studies, and for every export node its type, dataset,
// expressions and file name. Read-only: nothing is solved or saved.
//
//   comsol compile InspectMeshCell.java
//   comsol batch -inputfile InspectMeshCell.class
import com.comsol.model.*;
import com.comsol.model.util.*;

public class InspectMeshCell {
  public static void main(String[] args) throws java.io.IOException { run(); }

  public static Model run() throws java.io.IOException {
    // Fixed path: COMSOL's default security settings forbid reading system
    // properties from a batch class.
    String in = "/eos/user/a/akallits/p2_sim/comsol/MeshUnitCellIon_PressureStudy.mph";
    Model m = ModelUtil.load("Model", in);
    System.out.println("== parameters");
    for (String p : m.param().varnames())
      System.out.println("  " + p + " = " + m.param().get(p));
    System.out.println("== studies");
    for (String s : m.study().tags())
      System.out.println("  " + s + " : " + m.study(s).label());
    System.out.println("== physics");
    for (String ph : m.component("comp1").physics().tags())
      System.out.println("  " + ph + " : " + m.component("comp1").physics(ph).label());
    System.out.println("== datasets");
    for (String d : m.result().dataset().tags())
      System.out.println("  " + d + " : " + m.result().dataset(d).label());
    System.out.println("== exports");
    for (String e : m.result().export().tags()) {
      ExportFeature x = m.result().export(e);
      System.out.println("  " + e + " : " + x.label());
      for (String key : new String[] {"filename", "data", "expr", "descr", "unit",
                                      "struct", "header", "fullprec", "resolution", "location"}) {
        try {
          String[] v = x.getStringArray(key);
          System.out.println("      " + key + " = " + String.join(" | ", v));
        } catch (Exception ex) {
          try { System.out.println("      " + key + " = " + x.getString(key)); }
          catch (Exception ex2) { }
        }
      }
    }
    return m;
  }
}
