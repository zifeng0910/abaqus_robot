"""Build the single authorized 8.333 ms Wall-ON G6/L45 probe.

Only output cadence and step duration differ from the verified Wall-ON
reference deck; physical parameters and the mesh are copied unchanged.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003.inp"
JOB = "Wobble_F30_G6L45_WallOn_Free_0083_WobbleSurvival"
dst = ROOT / f"{JOB}.inp"
text = SRC.read_text(encoding="utf-8", errors="ignore")
text = re.sub(r"(?m)^\*Heading\s*$",
              f"*Heading\n** AUTHORIZED 8.333 ms WALL-ON WOBBLE-SURVIVAL PROBE: {JOB}; G6/L45; f=30 Hz",
              text, count=1)
text = re.sub(r"(?m)^\*Dynamic, Explicit\s*\r?\n,\s*0\.003\s*$",
              "*Dynamic, Explicit\n, 0.008333", text, count=1)
text = text.replace("** duration=0.003 s", "** duration=0.008333 s")
text = text.replace("** duration=0.03 s", "** duration=0.008333 s")
text = text.replace("duration=0.03 s", "duration=0.008333 s")
text = text.replace("duration=3 ms", "duration=8.333 ms")
text = text.replace("gradient 3 mT", "gradient 6 mT")
text = text.replace("gradient=3.000 mT=0.003000 T", "gradient=6.000 mT=0.006000 T")
text = text.replace("** gradient=6.000 mT=0.006000 T; server-side analytic field",
                    "** gradient=6.000 mT=0.006000 T; L=45 mm; server-side analytic field")
text = text.replace("** WALL/CEL FACTORIAL: WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003;",
                    f"** WALL/CEL FACTORIAL: {JOB};")
text = text.replace("** WALL-ON FREE-TRANSLATION PROBE: WobbleCal_F30_Cone30_B10_Grad6Forward_WallOn_Free_003;",
                    f"** WALL-ON FREE-TRANSLATION PROBE: {JOB};")
text = re.sub(r"time interval=1\.0e-4", "time interval=5.0e-5", text)
dst.write_text(text, encoding="utf-8")
print(f"WROTE {dst}")
