"""Build the single COM-fixed/free-rotation calibration deck from Probe080.

The robot reference node is relocated to the volume-weighted COM of the
existing Robot_SOLID tetrahedral mesh (including the assembly instance
translation). Translations are then constrained; rotations remain free.
"""
from pathlib import Path
import re
import numpy as np

ROOT = Path(__file__).resolve().parent
BASE = ROOT / 'CEL_R014_AnalyticRotField_Grad3mT_L45_D055_Bias40_Transverse90_Probe080_WallOff_CELOn.inp'
JOB = 'WobbleCal_COMFixed_F40_Cone30_B10'
OUT = ROOT / f'{JOB}.inp'


def robot_com_global(text: str):
    lines = text.splitlines()
    a = next(i for i, s in enumerate(lines) if s.startswith('*Part, name=Robot_SOLID'))
    en = next(i for i in range(a + 1, len(lines)) if lines[i].startswith('*Element'))
    ee = next(i for i in range(en + 1, len(lines)) if lines[i].startswith('*Elset') or lines[i].startswith('*Nset') or lines[i].startswith('*End Part'))
    nodes = {}
    for s in lines[a:en]:
        q = [x.strip() for x in s.strip().split(',')]
        if len(q) >= 4:
            try:
                nodes[int(q[0])] = np.array([float(q[1]), float(q[2]), float(q[3])])
            except ValueError:
                pass
    elems = []
    for s in lines[en + 1:ee]:
        q = [x.strip() for x in s.strip().split(',')]
        if len(q) >= 5:
            try:
                elems.append([int(x) for x in q[1:5]])
            except ValueError:
                pass
    vol = 0.0
    csum = np.zeros(3)
    for e in elems:
        x = [nodes[k] for k in e]
        ve = abs(np.linalg.det(np.column_stack((x[1] - x[0], x[2] - x[0], x[3] - x[0])))) / 6.0
        vol += ve
        csum += ve * sum(x) / 4.0
    c_local = csum / vol
    # The assembly instance translation is the only transform on Robot_SOLID.
    ia = text.index('*Instance, name=Robot_SOLID-1, part=Robot_SOLID')
    ie = text.index('*End Instance', ia)
    trans = None
    for s in text[ia:ie].splitlines()[1:]:
        q = [x.strip() for x in s.strip().split(',')]
        if len(q) >= 3:
            try:
                trans = np.array([float(q[0]), float(q[1]), float(q[2])])
                break
            except ValueError:
                pass
    if trans is None:
        raise RuntimeError('Robot instance translation not found')
    return c_local + trans, c_local, trans, len(nodes), len(elems), vol


def main():
    text = BASE.read_text(encoding='utf-8', errors='ignore')
    c_global, c_local, trans, n_nodes, n_elems, vol = robot_com_global(text)
    # Replace the assembly RP_ROBOT coordinate (the last *Node block before RP_PIPE).
    marker = '*Nset, nset=RP_PIPE'
    pre, post = text.split(marker, 1)
    # Replace the coordinate line belonging to assembly RP node 2 immediately
    # before RP_PIPE; do not rely on a fragile count of *Node blocks.
    tail_start = max(0, len(pre) - 2000)
    m = re.search(r'(?m)^([ \t]*2[ \t]*,[ \t]*)[-+0-9.eE]+[ \t]*,[ \t]*[-+0-9.eE]+[ \t]*,[ \t]*[-+0-9.eE]+[ \t]*$', pre[tail_start:])
    if not m:
        raise RuntimeError('assembly RP node 2 coordinate line not found')
    rpline_start = tail_start + m.start()
    rpline_end = tail_start + m.end()
    new_rp = '      2, {:.12f}, {:.12f}, {:.12f}'.format(*c_global)
    pre = pre[:rpline_start] + new_rp + pre[rpline_end:]
    text = pre + marker + post

    text = text.replace('CEL_R014_AnalyticRotField_Grad3mT_L45_D055_Bias40_Transverse90_Probe080_WallOff_CELOn', JOB)
    text = text.replace('** HighEnd83; FORWARD; Z90; lead5; XY(+2,-6); axis bias +40; 0.006 s',
                        '** COM-FIXED CALIBRATION: F40; cone30; B=10 mT; Wall-OFF; CEL-ON; translations fixed, rotations free')
    text = re.sub(r'(\*Dynamic, explicit\s*\n\s*,\s*)0\.0018', r'\g<1>0.075', text, count=1, flags=re.I)
    text = re.sub(r'(?i)2\.0e-6|2e-6', '1.0e-4', text)
    text = text.replace('** SOCKET MAGNETIC LOAD INTERFACE -- RP_ROBOT has all six DOF free.',
                        '** SOCKET MAGNETIC LOAD INTERFACE -- RP_ROBOT translations fixed at volume COM; UR1-UR3 free.')
    text = text.replace('** Six VUAMP amplitudes are updated by vuforc_socket_bridge.f.',
                        '** Six VUAMP amplitudes remain updated every Explicit increment by vuforc_socket_bridge.f.')
    text = text.replace('RP_PIPE, ENCASTRE\n', 'RP_PIPE, ENCASTRE\nRP_ROBOT, 1, 3\n', 1)
    text = text.replace('** RP = (4.90539503, -4.92436457, -6.52662897) mm',
                        '** COM-fixed RP_ROBOT = ({:.12f}, {:.12f}, {:.12f}) mm; volume-weighted mesh COM'.format(*c_global))
    text = text.replace('** Local centerline frame at RP (metadata only; no kinematic constraint)',
                        '** Local centerline frame at COM RP (metadata only); COM translations constrained in calibration')
    OUT.write_text(text, encoding='utf-8', newline='\n')
    print('created', OUT)
    print('robot_mesh_nodes', n_nodes, 'elements', n_elems, 'volume_mm3', vol)
    print('com_local_mm', c_local.tolist(), 'instance_translation_mm', trans.tolist(), 'com_global_mm', c_global.tolist())


if __name__ == '__main__':
    main()
