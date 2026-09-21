"""Report and render the terminated TAIL-gap diagnostic without extrapolating to 2T."""
from __future__ import annotations

import json
import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

import analyze_a14p5_tailgap_refine as analysis
import render_a14p5_tailgap_refine as render


OUT = analysis.OUT
JOB = analysis.FINE
END = 0.0125


def main():
    coarse = analysis.load_case(analysis.COARSE)
    fine = analysis.load_case(JOB)
    identity = fine["identity"]
    assert identity["status"] == "FAILED" and identity["dynamics_run_count"] == 1
    fields = render.visual_fields(fine)
    assert fields["time"][-1] >= END and fine["time"][-1] >= END
    closure0 = analysis.force_closure(coarse)["critical_8p0_12p5_ms"]
    closure1 = analysis.force_closure(fine)["critical_8p0_12p5_ms"]
    contact0, contact1 = coarse["contact"], fine["contact"]
    energy0, energy1 = analysis.energy_metrics(coarse), analysis.energy_metrics(fine)
    sta = (fine["case"] / (JOB + ".sta")).read_text(encoding="latin1")
    rows = re.findall(r"(?m)^\s*\d+\s+([0-9.E+-]+)\s+[0-9.E+-]+\s+\S+\s+([0-9.E+-]+)\s+\d+", sta)
    last_time = float(rows[-1][0])
    frame_end = float(fields["time"][-1])
    same_time = min(frame_end, fine["time"][-1])
    delta = float(np.interp(same_time, fine["time"], fine["axial"]) - fine["axial"][0])
    p = fields["pressure"]
    evf = fields["evf"]
    wet = (evf > 0.01) & np.isfinite(p)
    pclip = float(np.percentile(np.abs(p[wet]), 99))
    pclip = max(pclip, 1e-9)
    nodes, colors = render.robot_geometry(fine)
    radius = float(identity["lumen_radius_mm"])

    def make_gif(path, selector, fluid=False, zoom=False):
        indices = np.flatnonzero(selector)
        indices = np.unique(np.linspace(indices[0], indices[-1], min(100, len(indices))).astype(int))
        fig, ax = plt.subplots(figsize=(11, 3.8 if not zoom else 4.4))

        def draw(i):
            k = indices[i]
            t = float(fields["time"][k])
            ax.clear()
            if fluid:
                b = 0.0
                if zoom:
                    robot = render.robot_at(fine, nodes, t)
                    tail = robot[np.argmin((robot - fine["pipe"]) @ fine["c"])]
                    b = float((tail - fine["pipe"]) @ fine["b"])
                render.fluid_scatter(ax, fields, k, p, -pclip, pclip, b, zoom, zoom)
            render.draw_robot(ax, fine, nodes, colors, t)
            ax.axhline(radius, color="#5996a5")
            ax.axhline(-radius, color="#5996a5")
            if zoom:
                robot = render.robot_at(fine, nodes, t)
                tail_s = float(np.min((robot - fine["pipe"]) @ fine["c"]))
                ax.set_xlim(tail_s - 0.75, tail_s + 0.75)
            else:
                ax.set_xlim(-6, 6)
            ax.set_ylim(-0.8, 0.8)
            ax.set_aspect("equal")
            ax.set_xlabel("canonical +s (mm), left to right")
            ax.set_ylabel("n_routeA (mm)")
            ds = float(np.interp(t, fine["time"], fine["axial"]) - fine["axial"][0])
            title = "TAIL-gap pressure zoom" if zoom else ("TAIL-gap fluid" if fluid else JOB)
            ax.set_title(f"{title} | PARTIAL, TERMINATED at {last_time*1e3:.3f} ms / 16.667 ms\n"
                         f"t={t*1e3:.3f} ms | 120 Hz | B0=12 mT | G=2 mT | delta_s={ds:+.4f} mm")
            fig.tight_layout()

        FuncAnimation(fig, draw, frames=len(indices), interval=80).save(
            OUT / path, writer=PillowWriter(fps=12.5), dpi=100)
        plt.close(fig)
        return len(indices)

    gifs = {
        "fixed_side": "TRUECEL_A14P5_TAILGAP_REFINE_PARTIAL.gif",
        "fluid": "TRUECEL_A14P5_TAILGAP_REFINE_PARTIAL_FLUID.gif",
        "tail_zoom": "TRUECEL_A14P5_TAILGAP_REFINE_PARTIAL_TAIL_ZOOM.gif",
    }
    counts = {
        "fixed_side": make_gif(gifs["fixed_side"], fields["time"] <= frame_end),
        "fluid": make_gif(gifs["fluid"], fields["time"] <= frame_end, fluid=True),
        "tail_zoom": make_gif(gifs["tail_zoom"], (fields["time"] >= 0.008) & (fields["time"] <= END), fluid=True, zoom=True),
    }
    j0 = closure0["J_fluid_s_Ns_inferred_whole_minus_wall"]
    j1 = closure1["J_fluid_s_Ns_inferred_whole_minus_wall"]
    result = {
        "case": JOB,
        "status": "TERMINATED_BY_USER; PARTIAL; NOT_TWO_CYCLE_CONVERGENCE",
        "commanded_end_ms": 2000 / 120,
        "last_solver_time_ms": last_time * 1e3,
        "last_field_time_ms": frame_end * 1e3,
        "last_history_time_ms": float(fine["time"][-1] * 1e3),
        "available_delta_s_mm": delta,
        "common_impulse_window_ms": [8, 12.5],
        "coarse_J_fluid_s_Ns_inferred": j0,
        "refined_J_fluid_s_Ns_inferred": j1,
        "impulse_ratio_refined_over_coarse": j1 / j0,
        "refined_closure_error_fraction": closure1["relative_closure_error"],
        "contact_to_partial_end": contact1,
        "coarse_contact_full": contact0,
        "energy": {"coarse_ETOTAL": energy0["ETOTAL"], "refined_partial_ETOTAL": energy1["ETOTAL"]},
        "stable_dt_min_s": energy1["stable_timestep_s"]["min"],
        "general_contact_deep_penetration_warning": "InfoNodeDeepPenetFirst" in sta,
        "mesh_classification": "INDETERMINATE_INCOMPLETE_TWO_CYCLES",
        "recoil": "INDETERMINATE_INCOMPLETE_TWO_CYCLES",
        "contact": "WORSE_NUMERICAL_WARNING; episode comparison provisional",
        "numerics": "WORSE",
        "dynamics_run_count": 1,
        "mesh_type": identity["mesh_refinement_type"],
        "gifs": gifs,
        "gif_frame_counts": counts,
    }
    (OUT / (JOB + "_Partial_Analysis.json")).write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    lines = [
        f"# {JOB}: terminated partial diagnostic", "",
        "This is NOT a completed two-cycle convergence result. Abaqus was terminated at "
        f"{last_time*1e3:.3f} / {2000/120:.3f} ms on user request; last field frame {frame_end*1e3:.3f} ms.", "",
        "## Existing evidence", "",
        f"- Common 8.0-12.5 ms inferred robot-fluid axial impulse: coarse `{j0:+.7e}` N s; "
        f"refined `{j1:+.7e}` N s; ratio `{j1/j0:.4f}`.",
        f"- Refined momentum closure error in the common window: `{closure1['relative_closure_error']:.2%}`.",
        f"- Available displacement through {same_time*1e3:.3f} ms: `{delta:+.6f}` mm (not a two-cycle delta).",
        f"- Refined ETOTAL max absolute drift to stop: `{energy1['ETOTAL']['max_abs_drift_Nmm']:.6g}` N mm; "
        f"minimum reported stable increment `{energy1['stable_timestep_s']['min']:.3g}` s.",
        "- Abaqus warned that General Contact nodes penetrated tracked faces by over 50% of the "
        "typical 0.14599-mm element dimension (InfoNodeDeepPenetFirst).",
        "- Robot-fluid force is inferred as whole-robot General Contact minus direct robot-wall contact; "
        "it is not a native pair-isolated robot-fluid measurement.", "",
        "## Decision boundary", "",
        "Mesh classification and second-cycle recoil: INDETERMINATE because the run stopped before 2T. "
        "The 8.0-12.5 ms impulse is measured over a common window, but deep penetration and large "
        "energy drift make a physical mesh-convergence claim unsafe. NUMERICS: WORSE. "
        "This was a fixed-17,600-element local radial r-refinement, not added-cell h-refinement.", "",
        "All new GIFs are labeled PARTIAL/TERMINATED, use existing ODB frames only, and are not a completed scientific case.",
    ]
    (OUT / (JOB + "_Partial_Analysis.md")).write_text("\n".join(lines) + "\n", encoding="ascii")
    print(json.dumps({k: result[k] for k in ("last_field_time_ms", "refined_J_fluid_s_Ns_inferred", "impulse_ratio_refined_over_coarse", "mesh_classification", "gifs")}, indent=2))


if __name__ == "__main__":
    main()
