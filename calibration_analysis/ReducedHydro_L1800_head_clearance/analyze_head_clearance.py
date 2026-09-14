"""Run the validated exact-CAD analyzer against the HeadClearance job."""
from pathlib import Path
import importlib.util
import shutil


HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/"ReducedHydro_FreeCAD_L1800_validation"/"analyze_freecad_8p333.py"
JOB="Wobble_F30_G6L45_ReducedHydro_Zeta050_CAD_L1800_D0815_HeadClearance_0083"


def main():
    # The reused analyzer consumes this conventional surface-table name.
    shutil.copyfile(HERE/"headclearance_mesh060_surface_triangles.csv",
                    HERE/"freecad_robot_surface_triangles_exact.csv")
    spec=importlib.util.spec_from_file_location("validated_freecad_analyzer",SOURCE)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.HERE=HERE;module.JOB=JOB
    module.main()


if __name__=="__main__":main()
