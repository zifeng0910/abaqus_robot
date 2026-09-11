from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent
src = ROOT / 'WobbleCal_COMFixed_F30_Cone30_B10.inp'
job = sys.argv[1] if len(sys.argv) > 1 else 'WobbleCal_F30_Cone30_B10_WallOn_Free_003'
text = src.read_text(encoding='utf-8', errors='ignore')

text = re.sub(r'(?m)^\*Heading\s*$',
              '*Heading\n** WALL-ON FREE-TRANSLATION PROBE: %s; spin=30 Hz; duration=3 ms' % job,
              text, count=1)
text = text.replace('** WALL/CEL FACTORIAL: WobbleCal_COMFixed_F30_Cone30_B10; Wall=OFF;',
                    '** WALL/CEL FACTORIAL: %s; Wall=ON;' % job)
text = text.replace('** COM-FIXED CALIBRATION: F40; cone30; B=10 mT; Wall-OFF;',
                    '** FREE-TRANSLATION PROBE: F30; cone30; B=10 mT; Wall-ON;')
text = text.replace('** duration=0.03 s', '** duration=0.003 s')
text = re.sub(r'(?m)^,\s*0\.100\s*$', ', 0.003', text, count=1)

# Remove COM-fixed translational constraint; RP_PIPE remains encastre.
text = re.sub(r'(?m)^RP_ROBOT,\s*1,\s*3\s*\r?\n', '', text, count=1)
text = text.replace('** COM-FIXED CALIBRATION PARAMETERS:',
                    '** FREE-TRANSLATION PARAMETERS:')
text = text.replace('** SOCKET MAGNETIC LOAD INTERFACE -- RP_ROBOT translations fixed at volume COM; UR1-UR3 free.',
                    '** SOCKET MAGNETIC LOAD INTERFACE -- RP_ROBOT translations and UR1-UR3 free.')

# Add robot-to-wall contact while retaining robot/fluid and wall/fluid CEL
# coupling.  The helper wall surface is already present in this input.
needle = 'ROBOT_SOLID-1.ROBOT_SOLID_SURF, FLUID_CEL_SURF\n'
if 'ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF' not in text:
    text = text.replace('*Contact Inclusions\n' + needle,
                        '*Contact Inclusions\n'
                        'ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF\n'
                        + needle, 1)
    prop = 'ROBOT_SOLID-1.ROBOT_SOLID_SURF, FLUID_CEL_SURF, PROP_CEL_HARD\n'
    text = text.replace('*Contact Property Assignment\n' + prop,
                        '*Contact Property Assignment\n'
                        'ROBOT_SOLID-1.ROBOT_SOLID_SURF, Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF, PROP_CEL_HARD\n'
                        + prop, 1)

(ROOT / (job + '.inp')).write_text(text, encoding='utf-8')
print('WROTE', job + '.inp')
