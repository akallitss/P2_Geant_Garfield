#!/usr/bin/env python3
"""P2 final R&D phase -- simulation campaign + bulk tests (short deck).

Deliberately plain: white slides, a title, bullets, at most one figure.
Upload the .pptx to Google Drive and open it with Google Slides.
Regenerate with
    python docs/slides/make_rnd_slides.py
"""
import os

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE
from pptx.oxml import parse_xml
from pptx.util import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
MASK = os.path.join(ROOT, "Detector_Drawings", "Version_Sept26", "Insulation_masks",
                    "V3B_final", "QC", "P2_BASKET-Mask_M2_V3B_final_overview.png")
GAINFIG = os.path.join(HERE, "vmm_vs_dream_gain.jpg")   # from P2_Planning_2026_2027
OUT = os.path.join(HERE, "P2_RnD_final_phase_2026-10-07.pptx")

BLACK = RGBColor(0x20, 0x20, 0x20)
GREY = RGBColor(0x70, 0x70, 0x70)

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
BLANK = prs.slide_layouts[6]


def add_text(s, x, y, w, h, lines, size=20):
    """lines: '# heading', '- bullet', '  - sub-bullet', or plain text."""
    tf = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)).text_frame
    tf.word_wrap = True
    first = True
    for ln in lines:
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        if ln.startswith("# "):
            p.text, p.font.bold, p.font.size = ln[2:], True, Pt(size)
            p.space_before = Pt(10)
        elif ln.startswith("  - "):
            p.text, p.font.size, p.font.color.rgb = "– " + ln[4:], Pt(size - 4), GREY
            p.level = 1
        elif ln.startswith("- "):
            p.text, p.font.size = "• " + ln[2:], Pt(size)
        else:
            p.text, p.font.size = ln, Pt(size)
        if not ln.startswith("  - "):
            p.font.color.rgb = BLACK
        p.space_after = Pt(4)


def slide(title):
    s = prs.slides.add_slide(BLANK)
    tf = s.shapes.add_textbox(Inches(0.6), Inches(0.4), Inches(12.1), Inches(0.9)).text_frame
    p = tf.paragraphs[0]
    p.text, p.font.size, p.font.bold, p.font.color.rgb = title, Pt(30), True, BLACK
    return s


# 1 -- title
s = prs.slides.add_slide(BLANK)
add_text(s, 0.8, 2.6, 11.7, 1.0, ["P2 Micromegas — final R&D phase"], size=40)
add_text(s, 0.8, 3.6, 11.7, 1.0,
         ["Gain recovery tests, simulation campaign, new bulks"], size=24)
add_text(s, 0.8, 6.2, 11.7, 0.6,
         ["Alexandra Kallitsopoulou · CEA Saclay · 7 October 2026"], size=16)

# 2 -- why
s = slide("Why: VMM loses the low-gain part of the chamber")
add_text(s, 0.6, 1.4, 5.6, 5.8, [
    "- SPS July 2026, same pads, same beam:",
    "  - DREAM 95.6 %, VMM 85.3 %",
    "- The chamber gain varies by ×3.9 across the surface (non-planarity: drift gap and PCB bending)",
    "- Low-gain pads fall under the VMM threshold → lost",
    "# Plan: raise the gain and make it less sensitive to the gap",
    "- Larger isobutane fraction (Ar/iso 95/5 → 90/10)",
    "- Larger drift gap (4 → 5 mm)",
    "- Possibly other mixtures (Ne/iso 95/5, 85/15)",
    "- Measure the effect with the VMM electronics",
], size=18)
s.shapes.add_picture(GAINFIG, Inches(6.4), Inches(1.6), width=Inches(6.6))

# 3 -- simulation goals
s = slide("Simulation campaign: goals")
add_text(s, 0.6, 1.4, 12.1, 5.8, [
    "For each gas (Ar/iso 95/5, 90/10 · Ne/iso 95/5, 85/15) and drift gap (4, 5 mm):",
    "# 1. Gain and signal",
    "- Gain vs mesh voltage, and how much it changes with the local gap",
    "- Primary electrons per muon → charge on the pad → VMM efficiency (gain 3.0 / 4.5 mV/fC, 100 / 200 ns)",
    "- Check against the SPS data (DREAM and VMM)",
    "# 2. Background",
    "- Photons 40–100 keV: probability of ≥ 1 primary, Ne vs Ar",
    "- Event-by-event deposits: photons vs 100 MeV electrons → upper ADC cut",
    "- Weighted with Matthieu's background spectra (P2Sim) → rates per layer",
], size=19)

# 4 -- what we need
s = slide("What Geant4 + Garfield need")
add_text(s, 0.6, 1.4, 6.0, 5.8, [
    "# Geant4 (energy deposits)",
    "- Detector geometry from the design files ✓",
    "- Pillar map ✓",
    "- Gases with correct composition ✓",
    "- Particles: 150 GeV μ (SPS), 100 MeV e⁻, photons 10–150 keV ✓",
    "- Background spectra and directions (Matthieu) ✓",
], size=19)
add_text(s, 6.8, 1.4, 6.0, 5.8, [
    "# Garfield++ (charges and signals)",
    "- Gas tables (Magboltz) per mixture: drift velocity, diffusion, attachment, Townsend coefficient",
    "- Penning transfer rate per mixture (from measured gain curves)",
    "- Fields: drift and 150 µm amplification gap, mesh / drift voltages",
    "- Signal: pad weighting field → induced current",
    "- VMM shaper transfer function ✓ (Iakovidis et al.)",
    "- Measured gain curves (SPS HV scans) to calibrate",
], size=19)

# 5 -- steps
s = slide("Basic steps")
add_text(s, 0.6, 1.4, 6.4, 5.8, [
    "1.  Geant4: μ / e⁻ / γ through the P2 wedge → ionization clusters in the gas   ✓",
    "2.  Magboltz: gas tables for each mixture   → next",
    "3.  Drift to the mesh: diffusion, attachment, transparency",
    "4.  Avalanche: gain vs mesh voltage, calibrated on the SPS HV scans",
    "5.  Signal on the pad → VMM shaper → amplitude, time, threshold",
    "6.  Compare with SPS data, then predict 90/10, 5 mm, Ne",
    "7.  Background: photon runs × P2Sim spectra → e/γ separation, ADC cut",
], size=16)


def box(s, x, y, w, h, title, body, fill):
    shp = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    shp.fill.solid()
    shp.fill.fore_color.rgb = fill
    shp.line.color.rgb = GREY
    tf = shp.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text, p.font.bold, p.font.size, p.font.color.rgb = title, True, Pt(15), BLACK
    q = tf.add_paragraph()
    q.text, q.font.size, q.font.color.rgb = body, Pt(12), BLACK
    return shp


def arrow(s, x1, y1, x2, y2):
    c = s.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1), Inches(y1),
                               Inches(x2), Inches(y2))
    c.line.color.rgb = GREY
    c.line.width = Pt(2)
    c.line._get_or_add_ln().append(parse_xml(
        '<a:tailEnd type="triangle" xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"/>'))


WHITE = RGBColor(0xFF, 0xFF, 0xFF)
LIGHT = RGBColor(0xEE, 0xEE, 0xEE)
X, W = 7.3, 3.4
box(s, X, 1.5, W, 1.25, "Stage A · Geant4",
    "particle through the wedge → ionization clusters (position, time, primaries)  [step 1]",
    WHITE)
box(s, X, 3.3, W, 1.25, "Stage B · Digitizer",
    "drift → diffusion → attachment → transparency → gain → pad  [steps 3–4]", WHITE)
box(s, X, 5.1, W, 1.25, "Stage C · VMM3a",
    "shaper (Iakovidis) → threshold → amplitude, time  [step 5]", WHITE)
arrow(s, X + W / 2, 2.75, X + W / 2, 3.3)
arrow(s, X + W / 2, 4.55, X + W / 2, 5.1)
box(s, 11.05, 3.0, 2.0, 2.4, "Garfield++",
    "gas tables (Magboltz) · avalanche gain · pad signal  [steps 2, 4, 5]", LIGHT)
arrow(s, 11.05, 4.2, X + W, 3.93)
arrow(s, 11.05, 4.8, X + W, 5.5)

# 6 -- background (Matthieu)
s = slide("Background estimation (Matthieu, P2Sim)")
add_text(s, 0.6, 1.4, 12.1, 5.8, [
    "- Full P2 setup: photon rates and spectra at each Micromegas layer",
    "- ~65 % of the photons at 40–100 keV, mostly within 30° of the normal",
    "- Soft photons (< 20 keV) from the chamber wall: small flux, but half the detected rate in the first layer; stopped by its readout plane",
    "- Our campaign uses these spectra as input: response per photon from our detailed detector simulation × his rates",
    "- Geometries compared layer by layer; differences to align: mesh pitch, amplification gap, readout copper thickness",
], size=19)

# 7 -- bulk lab
s = slide("In parallel: new bulks in the lab")
add_text(s, 0.6, 1.4, 6.3, 5.8, [
    "# New insulation masks (V3B_final)",
    "- Pillars Ø 0.5 mm everywhere",
    "- 5 zones, pitch 2.0 mm (CERN reference) to 5.9 mm",
    "- Dead area 2.3–6.5 % per zone, 3.5 % overall",
    "- Test the efficiency per zone → choose the production pitch",
    "# PCBs with pillars",
    "- Mesh glued, no lamination: needs gluing tooling",
    "- Compare pre-stressed mesh (CERN supplier) vs standard",
], size=18)
s.shapes.add_picture(MASK, Inches(7.1), Inches(1.4), height=Inches(5.7))

prs.save(OUT)
print("wrote", OUT)
