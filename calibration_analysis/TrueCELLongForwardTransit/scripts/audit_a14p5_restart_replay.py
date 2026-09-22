"""Audit recovery reproducibility and actual prescribed-replay kinematics."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
JOB = "TRUECEL_A14P5_COARSETRAJ_CEL_HREFINE_REPLAY"
ORIGINAL = "TRUECEL_A14P5_DIAG_2CYCLES"
CASE = ROOT / "case" / JOB
PRIVATE = CASE / "private"
RESTART_LO = 0.0081
CRIT_LO = 1.0 / 120.0
CRIT_HI = 0.010682020978413481


def split_index(t):
    jumps = np.flatnonzero(np.diff(t) < 0)
    if len(jumps) != 1:
        raise RuntimeError(f"Expected one recovery time reset, found {len(jumps)}")
    return int(jumps[0] + 1)


def vectors(z, stem):
    return z[f"{stem}1"][:, 0], np.column_stack([z[f"{stem}{i}"][:, 1] for i in (1, 2, 3)])


def stats(t, error):
    magnitude = np.linalg.norm(error, axis=1) if error.ndim == 2 else np.abs(error)
    i = int(np.argmax(magnitude))
    return {"rms": float(np.sqrt(np.mean(magnitude**2))), "p99": float(np.percentile(magnitude, 99)),
            "max": float(magnitude[i]), "max_time_ms": float(t[i] * 1e3)}


def compare_segments(t, y, q, lo, hi):
    old_t, new_t = t[:q], t[q:]
    mask = (new_t >= lo) & (new_t <= hi)
    reference = np.column_stack([np.interp(new_t[mask], old_t, y[:q, i]) for i in range(y.shape[1])])
    return new_t[mask], y[q:][mask] - reference


def integrate(t, y, lo, hi):
    inside = (t > lo) & (t < hi)
    x = np.r_[lo, t[inside], hi]
    return float(np.trapz(np.interp(x, t, y), x))


def contact_vector(z, suffix, prefix="CFT"):
    arrays = []
    t = None
    for axis in (1, 2, 3):
        key = next(k for k in z.files if f"|{prefix}{axis} " in k and k.endswith(suffix))
        t = z[key][:, 0]
        arrays.append(z[key][:, 1])
    return t, np.column_stack(arrays)


def authoritative(t, y, q):
    keep_old = t[:q] < t[q]
    return np.r_[t[:q][keep_old], t[q:]], np.vstack((y[:q][keep_old], y[q:]))


def main():
    identity = json.loads((CASE / "case_identity.json").read_text())
    s = np.asarray(identity["canonical_plus_s_axis_aba"])
    b = np.asarray(identity["b_routeA_aba"])
    rp = np.load(PRIVATE / "rp_history_private.npz")
    rt, u_all = vectors(rp, "U")
    q = split_index(rt)
    overlap_hi = float(rt[q - 1])
    overlap = {}
    authoritative_rp = {}
    for stem in ("U", "UR", "V", "VR"):
        t, y = vectors(rp, stem)
        if split_index(t) != q:
            raise RuntimeError(f"Inconsistent recovery split in {stem}")
        x, error = compare_segments(t, y, q, float(t[q]), overlap_hi)
        overlap[stem] = stats(x, error)
        at, ay = authoritative(t, y, q)
        for axis in range(3):
            authoritative_rp[f"{stem}{axis+1}"] = np.column_stack((at, ay[:, axis]))
    np.savez_compressed(PRIVATE / "rp_history_authoritative_private.npz", **authoritative_rp)

    contact = np.load(PRIVATE / "force_flux_contact_history_private.npz")
    whole_name = "on surface ASSEMBLY_ROBOT_SOLID-1_ROBOT_SOLID_SURF"
    wall_name = whole_name + "/ASSEMBLY_PIPE_WALL_HELPER-1_PIPE_WALL_HELPER_SURF"
    ct, whole = contact_vector(contact, whole_name)
    _, wall = contact_vector(contact, wall_name)
    cq = split_index(ct)
    contact_stats = {}
    for name, values in (("whole_robot", whole), ("direct_robot_wall", wall)):
        x, error = compare_segments(ct, values, cq, float(ct[cq]), float(ct[cq-1]))
        result = stats(x, error)
        old_scale = np.percentile(np.linalg.norm(values[:cq][(ct[:cq] >= ct[cq]) & (ct[:cq] <= ct[cq-1])], axis=1), 99)
        result["p99_relative_to_original_p99"] = result["p99"] / max(float(old_scale), 1e-15)
        old_j = integrate(ct[:cq], values[:cq] @ s, float(ct[cq]), float(ct[cq-1]))
        new_j = integrate(ct[cq:], values[cq:] @ s, float(ct[cq]), float(ct[cq-1]))
        result.update(axial_impulse_original_Ns=old_j, axial_impulse_recovered_Ns=new_j,
                      axial_impulse_relative_difference=abs(new_j-old_j)/max(abs(old_j), 1e-15))
        contact_stats[name] = result
    old_fluid = (whole[:cq] - wall[:cq]) @ s
    new_fluid = (whole[cq:] - wall[cq:]) @ s
    old_j = integrate(ct[:cq], old_fluid, float(ct[cq]), float(ct[cq-1]))
    new_j = integrate(ct[cq:], new_fluid, float(ct[cq]), float(ct[cq-1]))
    contact_stats["inferred_robot_CEL_axial_impulse"] = {
        "original_Ns": old_j, "recovered_Ns": new_j,
        "relative_difference": abs(new_j-old_j)/max(abs(old_j), 1e-15),
        "note": "Integrated on each native time grid; pointwise subtraction is cancellation-sensitive."
    }

    magnetic = np.loadtxt(CASE / "magnetic_increment.csv", delimiter=",", skiprows=1)
    expected_phase = np.mod(43200.0 * magnetic[:, 0], 360.0)
    phase_error = np.abs(((magnetic[:, 2] - expected_phase + 180) % 360) - 180)
    phase = {
        "first_time_ms": float(magnetic[0, 0] * 1e3), "last_time_ms": float(magnetic[-1, 0] * 1e3),
        "first_phase_deg": float(magnetic[0, 2]), "expected_first_phase_deg": float(expected_phase[0]),
        "max_formula_phase_error_deg": float(phase_error.max()),
        "total_time_not_reset": bool(magnetic[0, 0] > 0.008 and magnetic[-1, 0] >= 0.0108 - 1e-12),
    }
    fields = np.load(PRIVATE / "truecel_field_private.npz")
    field_time = fields["time"]
    field_note = ("Recovered ODB contains one monotonic field-frame sequence; Abaqus replaced/truncated the old overlap "
                  "frames, so two independent CEL field histories are not available for overlap comparison.")

    passed = (overlap["U"]["max"] < 1e-6 and overlap["UR"]["max"] < 1e-6 and
              overlap["V"]["p99"] < 0.05 and overlap["VR"]["p99"] < 0.25 and
              all(contact_stats[k]["p99_relative_to_original_p99"] < 0.10 and
                  contact_stats[k]["axial_impulse_relative_difference"] < 0.01
                  for k in ("whole_robot", "direct_robot_wall")) and
              contact_stats["inferred_robot_CEL_axial_impulse"]["relative_difference"] < 0.01 and
              phase["total_time_not_reset"] and phase["max_formula_phase_error_deg"] < 1e-6)
    classification = "RESTART_OVERLAP_REPRODUCES_ORIGINAL" if passed else "RESTART_OVERLAP_MISMATCH"
    restart_report = {
        "classification": classification, "overlap_ms": [float(rt[q]*1e3), overlap_hi*1e3],
        "kinematic_vector_norm_differences": overlap, "contact_resultant_differences": contact_stats,
        "absolute_time_and_VUAMP_phase": phase,
        "CEL_field_overlap": {"available": False, "field_time_monotonic": bool(np.all(np.diff(field_time) >= 0)), "note": field_note},
        "authoritative_handoff": "original segment for t < recovered first timestamp; recovered segment from that timestamp onward; no averaging or duplicate timestamps",
        "thresholds": {"U_max_mm": 1e-6, "UR_max_rad": 1e-6, "V_p99_mm_s": .05, "VR_p99_rad_s": .25,
                       "contact_p99_relative": .10, "contact_axial_impulse_relative": .01},
    }
    (ROOT / "A14P5_RESTART_REPRODUCIBILITY_AUDIT.json").write_text(json.dumps(restart_report, indent=2)+"\n")
    md = ["# A14P5 restart reproducibility audit", "", f"**{classification}**", "",
          f"Overlap: {rt[q]*1e3:.6f}-{overlap_hi*1e3:.6f} ms.", "",
          f"Kinematic vector-norm differences: `{json.dumps(overlap)}`", "",
          f"Contact resultant differences: `{json.dumps(contact_stats)}`", "",
          f"Absolute-time / VUAMP phase: `{json.dumps(phase)}`", "", field_note, "",
          "Authoritative handoff: original segment before the first recovered timestamp, recovered segment from that timestamp onward; no averaging and no duplicate timestamps."]
    (ROOT / "A14P5_RESTART_REPRODUCIBILITY_AUDIT.md").write_text("\n".join(md)+"\n")
    if not passed:
        print(json.dumps(restart_report, indent=2)); return

    auth = np.load(PRIVATE / "rp_history_authoritative_private.npz")
    source = np.load(ROOT / "case" / ORIGINAL / "private" / "rp_history_private.npz")
    t = auth["U1"][:, 0]
    mask = (t >= CRIT_LO) & (t <= CRIT_HI)
    tx = t[mask]
    def matrix(z, stem, target):
        base_t = z[f"{stem}1"][:, 0]
        return np.column_stack([np.interp(target, base_t, z[f"{stem}{i}"][:, 1]) for i in (1,2,3)])
    au, ar, av, aw = [matrix(auth, stem, tx) for stem in ("U","UR","V","VR")]
    su, sr, sv, sw = [matrix(source, stem, tx) for stem in ("U","UR","V","VR")]
    trans = np.linalg.norm(au-su, axis=1)
    orient = np.rad2deg((Rotation.from_rotvec(au*0+ar).inv()*Rotation.from_rotvec(sr)).magnitude())
    axial_v = np.abs((av-sv) @ s)
    rocking_w = np.abs((aw-sw) @ b)
    replay = {"translation_max_error_mm": float(trans.max()), "orientation_max_error_deg": float(orient.max()),
              "v_s_error_mm_s": stats(tx, axial_v), "rocking_axis_VR_error_rad_s": stats(tx, rocking_w)}
    # Force-event error samples.
    act, awhole = authoritative(ct, whole, cq); _, awall = authoritative(ct, wall, cq)
    fluid = (awhole-awall) @ s
    fm = (act >= CRIT_LO) & (act <= CRIT_HI)
    fi = np.flatnonzero(fm)
    peak = fi[np.argmin(fluid[fm])]
    deriv = np.abs(np.diff(fluid) / np.diff(act)); candidates = fi[1:][np.argsort(deriv[fi[:-1]])[-3:]]
    events = [peak] + list(candidates)
    replay["force_event_errors"] = [{"time_ms": float(act[i]*1e3), "event": "peak_negative_robot_CEL" if j==0 else "major_force_change",
        "v_s_error_mm_s": float(np.interp(act[i],tx,axial_v)), "rocking_axis_VR_error_rad_s": float(np.interp(act[i],tx,rocking_w)),
        "F_robotCEL_axial_N": float(fluid[i])} for j,i in enumerate(events)]
    replay["passed"] = bool(replay["translation_max_error_mm"] < .001 and replay["orientation_max_error_deg"] < .02)
    (ROOT / "A14P5_REFINED_REPLAY_VALIDATION.json").write_text(json.dumps(replay, indent=2)+"\n")
    print(json.dumps({"restart": restart_report, "replay": replay}, indent=2))


if __name__ == "__main__":
    main()
