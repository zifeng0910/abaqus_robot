"""Generate SCREENING ONLY straight-body-length variants with the frozen HEAD."""

from pathlib import Path
import csv
import hashlib
import json
import math

import FreeCAD as App
import Part


HERE = Path(__file__).resolve().parent
LENGTHS = (2.00, 2.10, 2.20, 2.30, 2.40)
RADIUS = 0.4075
HEAD_END = 0.26149477
DENSITY = 7.80906654321e-9
POLES = ((0.0, 0.0), (0.0, 0.09664525), (0.005, 0.18497527),
         (0.17271875, RADIUS), (0.25, RADIUS), (HEAD_END, RADIUS))


def shape_for(length):
    curve = Part.BezierCurve()
    curve.setPoles([App.Vector(x, r, 0.0) for x, r in POLES])
    edges = [curve.toShape(),
             Part.makeLine(App.Vector(HEAD_END, RADIUS, 0), App.Vector(length, RADIUS, 0)),
             Part.makeLine(App.Vector(length, RADIUS, 0), App.Vector(length, 0, 0)),
             Part.makeLine(App.Vector(length, 0, 0), App.Vector(0, 0, 0))]
    return Part.Face(Part.Wire(edges)).revolve(App.Vector(), App.Vector(1, 0, 0), 360).removeSplitter()


def profile_rows(length):
    rows = []
    degree = len(POLES) - 1
    for i in range(501):
        t = i / 500.0
        weights = [math.comb(degree, k) * (1-t)**(degree-k) * t**k for k in range(degree+1)]
        rows.append(("HEAD", t, sum(w*p[0] for w, p in zip(weights, POLES)),
                     sum(w*p[1] for w, p in zip(weights, POLES))))
    rows.extend(("BODY", "", x, RADIUS) for x in (HEAD_END, length))
    rows.extend(("TAIL", "", length, r) for r in (RADIUS, 0.0))
    return rows


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    manifest = []
    for length in LENGTHS:
        code = "%04d" % round(length * 1000)
        stem = "Robot_SCREENING_L%s_D0815" % code
        shape = shape_for(length)
        if not shape.isValid() or not shape.isClosed() or len(shape.Solids) != 1:
            raise RuntimeError("Invalid screening shape %s" % stem)
        box = shape.BoundBox
        if max(abs(box.XLength-length), abs(box.YLength-.815), abs(box.ZLength-.815)) > 1e-6:
            raise RuntimeError("Dimension failure %s" % stem)
        step = HERE / (stem + ".step")
        doc = App.newDocument(stem)
        obj = doc.addObject("PartDesign::Feature", "Robot")
        obj.Label = "SCREENING ONLY L%.2f D0.815 robot" % length
        obj.Shape = shape
        Part.export([obj], str(step))
        App.closeDocument(doc.Name)
        props = {
            "designation": "SCREENING ONLY - NOT AN AUTHORITATIVE DESIGN",
            "length_mm": length, "diameter_mm": .815,
            "head_definition": "frozen degree-5 Bezier HEAD from L2300 authoritative CAD",
            "volume_mm3": shape.Volume, "mass_mg": shape.Volume * DENSITY * 1e9,
            "center_of_mass_mm": [shape.CenterOfMass.x, shape.CenterOfMass.y, shape.CenterOfMass.z],
            "principal_inertia_tonne_mm2": [x*DENSITY for x in shape.PrincipalProperties["Moments"]],
            "step_sha256": hashlib.sha256(step.read_bytes()).hexdigest(),
        }
        (HERE / (stem + "_geometry.json")).write_text(json.dumps(props, indent=2) + "\n")
        with (HERE / (stem + "_profile.csv")).open("w", newline="") as handle:
            writer = csv.writer(handle); writer.writerow(("region", "parameter", "x_mm", "radius_mm"))
            writer.writerows(profile_rows(length))
        manifest.append(props)
    (HERE / "screening_geometry_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

