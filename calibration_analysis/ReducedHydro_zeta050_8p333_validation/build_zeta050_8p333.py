"""Build the frozen zeta=0.50 8.333 ms validation deck; never launches Abaqus."""
from pathlib import Path
import hashlib
import json


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORKDIR = REPO.parent
SOURCE_DIR = REPO / 'calibration_analysis' / 'ReducedHydro_hidden_impact_audit' / 'normal_contact_damping_probe'
SOURCE_JOB = 'Wobble_F30_G6L45_ReducedHydroFixed_NormalDamp050_ContactAudit_0013'
JOB = 'Wobble_F30_G6L45_ReducedHydro_Zeta050_WallOn_Free_0083'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    source_path = SOURCE_DIR / f'{SOURCE_JOB}.inp'
    source = source_path.read_text()
    assert source.count('1.0e-7, 0.001300') == 1
    assert source.count('*Output, field, time interval=5.0e-7, time marks=NO') == 1
    candidate = source.replace('1.0e-7, 0.001300', '1.0e-7, 0.008333', 1)
    candidate = candidate.replace(
        '*Output, field, time interval=5.0e-7, time marks=NO',
        '*Output, field, time interval=2.5e-5, time marks=NO', 1)
    candidate = candidate.replace(
        'save contact and energy fields every 1 us.',
        'save visualization/contact fields every 25 us; RP and contact-resultant history remains every increment.')

    frozen_tokens = [
        '*Contact Damping, definition=CRITICAL DAMPING FRACTION, tangent fraction=0.0\n0.50',
        '0.03,', '*Dynamic, Explicit, DIRECT\n1.0e-7, 0.008333',
        '*Output, history, frequency=1',
        '*Contact Output, surface=ROBOT_SOLID-1.ROBOT_SOLID_SURF\nCFN, CFS, CFT',
        '*Contact Output, surface=Pipe_WALL_HELPER-1.PIPE_WALL_HELPER_SURF\nCFN, CFS, CFT',
        'phase=248.0 deg', 'gradient=6.000 mT=0.006000 T',
    ]
    for token in frozen_tokens:
        assert token in candidate, token
    for variable in ('U1', 'U2', 'U3', 'UR1', 'UR2', 'UR3',
                     'V1', 'V2', 'V3', 'VR1', 'VR2', 'VR3'):
        assert f'name=RP_{variable}' in candidate
    bridge = REPO / 'calibration_analysis' / 'Wobble30Hz_reduced_hydrodynamics' / 'vuforc_socket_bridge_reduced_hydro.f'
    bridge_text = bridge.read_text()
    for token in ('cpar /4.0D-9/', 'cperp /1.2D-8/', 'kspin /1.0D-9/', 'kwob /3.0D-9/'):
        assert token in bridge_text, token

    local = HERE / f'{JOB}.inp'
    work = WORKDIR / f'{JOB}.inp'
    local.write_text(candidate)
    work.write_text(candidate)
    identity = {
        'job': JOB,
        'source_job': SOURCE_JOB,
        'source_sha256': sha(source_path),
        'candidate_sha256': sha(local),
        'duration_s': 0.008333,
        'direct_dt_s': 1e-7,
        'field_interval_s': 2.5e-5,
        'history_interval': 'every explicit increment',
        'allowed_nonphysics_changes': ['duration 0.001300 -> 0.008333 s',
                                       'field interval 5e-7 -> 2.5e-5 s',
                                       'output comment only'],
        'frozen': {
            'critical_damping_fraction': 0.50,
            'tangent_fraction': 0.0,
            'friction': 0.03,
            'B0_T': 0.010,
            'frequency_Hz': 30.0,
            'cone_deg': 30.0,
            'bias_deg': 40.0,
            'phase_deg': 248.0,
            'sense': 1,
            'gradient_T': 0.006,
            'gradient_length_mm': 45.0,
            'Cparallel_Ns_per_mm': 4e-9,
            'Cperp_Ns_per_mm': 1.2e-8,
            'Kspin_Nmm_s': 1e-9,
            'Kwobble_Nmm_s': 3e-9,
        },
    }
    (HERE / 'zeta050_8p333_input_identity.json').write_text(json.dumps(identity, indent=2))
    print(json.dumps(identity, indent=2))


if __name__ == '__main__':
    main()
