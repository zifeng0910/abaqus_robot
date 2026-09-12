"""Generate the fixed-camera animation for the L=1.800 mm candidate."""
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "ReducedHydro_zeta050_8p333_validation"
sys.path.insert(0, str(BASE))
import make_zeta050_gifs as animation


def main():
    animation.HERE = HERE
    animation.JOB = "Wobble_F30_G6L45_ReducedHydro_Zeta050_L1800_D0815_WallOn_Free_0083"
    animation.PREFIX = "L1800_8p333"
    animation.MAIN_GIF = "Wobble_F30_ReducedHydro_Zeta050_L1800_D0815_WallOn_8p333.gif"
    animation.MAKE_COMPARISON = False
    animation.main()


if __name__ == "__main__":
    main()
