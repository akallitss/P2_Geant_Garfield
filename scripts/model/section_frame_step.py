"""Section the P2 drift-frame STEP and print its profile — the source of the
frame constants in include/P2Wedge.hh (docs/P2_GEOMETRY.md §5).

Needs FreeCAD (OpenCASCADE); run with its interpreter, not plain python:

    freecad.cmd scripts/model/section_frame_step.py [path/to/frame.stp]

(The snap build cannot see /tmp; keep the script and the STEP under $HOME.)

Prints the z of every plane face, then the material intervals along radial
lines (r) and along arcs (phi) at a ladder of heights above the board face.
"""
import math
import sys

import Part
from FreeCAD import Vector as V

DEFAULT = ("/local/home/ak271430/Documents/PostDocSaclay/Detector_Drawings/"
           "DRIFT_P2_261124/DRIFT_P2_261124/FRAME_4mm/P2_Frame_V2.0.stp")
args = [a for a in sys.argv[1:] if a.endswith((".stp", ".step"))]
shape = Part.read(args[0] if args else DEFAULT)

print("bbox", shape.BoundBox, " solids", len(shape.Solids))
zs = {}
for f in shape.Faces:
    try:
        n = f.normalAt(0, 0)
    except Exception:
        continue
    if abs(abs(n.z) - 1) < 1e-6:
        z = round(f.BoundBox.ZMin, 4)
        zs[z] = zs.get(z, 0.0) + f.Area
for z in sorted(zs):
    print(f"plane z = {z:7.3f} mm   area {zs[z]:9.1f} mm2")

Z = [0.1, 1.0, 2.0, 2.45, 3.0, 3.95, 4.05, 5.0, 7.9]
for phi in (7.0, 30.0, 53.0):
    c, s = math.cos(math.radians(phi)), math.sin(math.radians(phi))
    print(f"\nradial line phi = {phi:g} deg: material r-intervals [mm]")
    for z in Z:
        com = shape.common(Part.makeLine(V(0, 0, z), V(800 * c, 800 * s, z)))
        iv = sorted(tuple(sorted(round(math.hypot(v.X, v.Y), 2)
                                 for v in (e.Vertexes[0], e.Vertexes[-1])))
                    for e in com.Edges)
        print(f"  z = {z:4.2f}  {iv}")
for R in (200.0, 550.0):
    print(f"\narc r = {R:g} mm: material phi-intervals [deg]")
    for z in Z:
        com = shape.common(Part.makeCircle(R, V(0, 0, z), V(0, 0, 1), -20, 80))
        iv = sorted(tuple(sorted(round(math.degrees(math.atan2(v.Y, v.X)), 3)
                                 for v in (e.Vertexes[0], e.Vertexes[-1])))
                    for e in com.Edges)
        print(f"  z = {z:4.2f}  {iv}")
