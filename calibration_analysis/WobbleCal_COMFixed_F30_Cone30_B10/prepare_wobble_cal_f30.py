"""Create the single-variable 30 Hz COM-fixed calibration deck."""
from pathlib import Path

root = Path(__file__).resolve().parent
src = root / "WobbleCal_COMFixed_F40_Cone30_B10.inp"
dst = root / "WobbleCal_COMFixed_F30_Cone30_B10.inp"
text = src.read_text(encoding="latin1")
text = text.replace("WobbleCal_COMFixed_F40_Cone30_B10", "WobbleCal_COMFixed_F30_Cone30_B10")
text = text.replace("phase=248 deg; gradient=3.000 mT=0.003000 T", "phase=248 deg; gradient=3.000 mT=0.003000 T; spin=30 Hz")
text = text.replace("\n, 0.075\n", "\n, 0.100\n")
dst.write_text(text, encoding="latin1")
print("WROTE", dst)
