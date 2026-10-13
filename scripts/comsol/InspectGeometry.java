// InspectGeometry.java -- print the parameters and every geometry feature
// (with all its properties) of the colleague's mesh-cell model, to trace how
// meshHeight, rWire and the weave define the wire positions.
//   comsol compile InspectGeometry.java
//   comsol batch -nosave -inputfile InspectGeometry.class
import com.comsol.model.*;
import com.comsol.model.util.*;

public class InspectGeometry {
  static final String IN =
      "/eos/user/a/akallits/p2_sim/comsol/MeshUnitCellIon_PressureStudy.mph";

  public static void main(String[] args) throws java.io.IOException { run(); }

  public static Model run() throws java.io.IOException {
    Model m = ModelUtil.load("Model", IN);
    System.out.println("== parameters");
    for (String p : m.param().varnames())
      System.out.println("  " + p + " = " + m.param().get(p) + "   // " + m.param().descr(p));
    for (String g : m.component("comp1").geom().tags()) {
      GeomSequence gs = m.component("comp1").geom(g);
      System.out.println("== geometry " + g + " (length unit " + gs.lengthUnit() + ")");
      for (String f : gs.feature().tags()) {
        GeomFeature ft = gs.feature(f);
        System.out.println("  " + f + " : " + ft.label() + " [" + ft.getType() + "]");
        for (String key : ft.properties()) {
          String val;
          try { val = String.join(" | ", ft.getStringArray(key)); }
          catch (Exception e1) {
            try { val = ft.getString(key); } catch (Exception e2) { continue; }
          }
          if (val.length() > 0 && val.length() < 300)
            System.out.println("      " + key + " = " + val);
        }
      }
    }
    return null;
  }
}
