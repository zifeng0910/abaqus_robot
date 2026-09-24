"""Prepare the sole fresh-start F80 TRUE-CEL candidate from the solved F100 deck."""
from __future__ import annotations

import hashlib, json, shutil
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OLD='TRUECEL_B0P11_G2P20_A14P5_F100_FAST'
NEW='TRUECEL_B0P11_G2P20_A14P5_F80_FAST'
SRC=ROOT/'case'/OLD; DST=ROOT/'case'/NEW

def once(s,a,b):
    if s.count(a)!=1: raise RuntimeError(f'Expected one {a!r}; got {s.count(a)}')
    return s.replace(a,b,1)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest().upper()

def main():
    if DST.exists():raise RuntimeError(f'Candidate exists: {DST}')
    ident=json.loads((SRC/'case_identity.json').read_text())
    for k,v in {'frequency_Hz':100.0,'B0_mT':11.0,'gradient_mT':2.2,
                'rocking_main_amplitude_deg':14.5,'rocking_cross_amplitude_deg':2.5,
                'fluid_EOS_c0_mm_s':100000.0,'Eulerian_dimensions':[44,20,20],
                'explicit_stable_time_scale_factor':0.4}.items():
        if ident.get(k)!=v:raise RuntimeError(f'{k} mismatch')
    inp=(SRC/f'{OLD}.inp').read_text(encoding='latin1').replace(OLD,NEW)
    inp=once(inp,', 0.020000000000',', 0.037500000000')
    f90=(SRC/'vuamp_precomputed_truecel.f90').read_text(encoding='latin1').replace(OLD,NEW)
    f90=once(f90,'phase_deg=modulo(36000.0d0*t,360.0d0)',
             'phase_deg=modulo(28800.0d0*t,360.0d0)')
    f90=f90.replace('monitor_next_boundary=1.0d-2','monitor_next_boundary=1.25d-2')
    f90=f90.replace('monitor_next_boundary+1.0d-2','monitor_next_boundary+1.25d-2')
    f90=f90.replace('5.0d-3','6.25d-3')
    f90=f90.replace('F100_RECOIL','F80_RECOIL')
    f90=f90.replace('f100_event.txt','f80_event.txt')
    f90=f90.replace('magnetic_increment_g2p20_f100.csv','magnetic_increment_g2p20_f80.csv')
    if f90.count('grad=0.002200d0*grad')!=1 or f90.count('b=0.011d0*b')!=1:
        raise RuntimeError('Magnetic amplitude changed')
    DST.mkdir()
    table='magnetic_field_gradient_table_B0P11_A14P5.dat'
    shutil.copy2(SRC/table,DST/table);shutil.copy2(SRC/(table+'.json'),DST/(table+'.json'))
    inp_path=DST/f'{NEW}.inp';f90_path=DST/'vuamp_precomputed_truecel.f90'
    inp_path.write_text(inp,encoding='latin1');f90_path.write_text(f90,encoding='latin1')
    for k in ('wallclock_s','cpus','socket_calls','failure'):ident.pop(k,None)
    ident.update({'case_id':NEW,'status':'PREPARED','classification':'F80_THREE_CYCLE_GATE',
                  'physical_parent':OLD,'frequency_Hz':80.0,'duration_s':.0375,
                  'dynamics_run_count':0,'single_physics_change':'frequency 100 -> 80 Hz',
                  'only_new_physics_change':'frequency 100 -> 80 Hz',
                  'period_s':.0125,'stage1_duration_s':.0375,
                  'input_sha256':sha(inp_path),'fortran_sha256':sha(f90_path),
                  'magnetic_table_sha256':sha(DST/table)})
    ident['frozen_for_candidate']=[v.replace('f=100Hz','f=80Hz') for v in ident['frozen_for_candidate']]
    ident['online_fail_fast']['position']='backtrack >0.020 mm and remains >0.010 mm below maximum for >=6.25 ms'
    (DST/'case_identity.json').write_text(json.dumps(ident,indent=2)+'\n')
    audit={'case':NEW,'parent':OLD,'fresh_start':True,'single_physics_change':'frequency 100 -> 80 Hz',
           'B0_mT':11,'G_mT':2.2,'frequency_Hz':80,'A_main_deg':14.5,'A_cross_deg':2.5,
           'duration_s':.0375,'mesh':[44,20,20],'c0_mm_s':100000,
           'input_sha256':sha(inp_path),'fortran_sha256':sha(f90_path),'table_sha256':sha(DST/table)}
    (ROOT/f'{NEW}_Setup_Audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))

if __name__=='__main__':main()
