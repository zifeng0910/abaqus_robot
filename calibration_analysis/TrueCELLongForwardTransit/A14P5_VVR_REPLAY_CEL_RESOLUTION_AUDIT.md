# A14P5 V/VR replay CEL resolution audit

Old UR replay is invalid (`FIXED_TRAJECTORY_REPLAY_INVALID`) and excluded. The global translational V(t) + spatial/global VR(t) driver was independently validated.

Case B coarse replay and Case C refined replay both passed the rotation-matrix/geodesic kinematic gates. A–B is replay-method bias only: S_replay = 1.40018. B–C is the authoritative CEL-resolution comparison: S_CEL = 0.755597.

Critical-window inferred impulses (8.333333–10.682021 ms): A=1.493524511e-07 N s, B=-5.976854493e-08 N s, C=-1.460761760e-08 N s. All five requested windows are in `A14P5_VVR_REPLAY_CEL_RESOLUTION_WINDOWS.csv`. Fixed-time B–C field RMS (native element ordering; pressure then EVF at 8.5, 9.0, 9.5, 10.0, 10.3, 10.682021 ms): {'p': [0.02126608043909073, 0.017927709966897964, 0.01379033550620079, 0.013237318955361843, 0.01277141459286213, 0.011007717810571194], 'e': [0.668285071849823, 0.6670253276824951, 0.6665053367614746, 0.6672095656394958, 0.6684111952781677, 0.6692406535148621]}. The synchronized GIF uses shared physical coordinates and pressure color limits.

Classification: `VVR_COARSE_AND_REFINED_REPLAY_VALID`.
CEL classification: `ROBOT_CEL_LOAD_STRONGLY_CEL_RESOLUTION_SENSITIVE`.

formal_mesh_convergence = NO; physical_fluid_load_validated = NO; CONTROL_VOLUME_UNAVAILABLE. Prescribed-motion constraint work means free-body ETOTAL propulsion criteria are not applied. Force labels are inferred robot-CEL = whole contact minus direct wall contact, not direct surface traction.
