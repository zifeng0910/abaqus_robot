# Human-in-the-loop coarse motion-mode screen

This directory contains a 32-case, one-factor-plus-directed-combo Abaqus/Explicit
screen. It deliberately produces no overall score and selects no winner. The common
baseline is `GEO_230`; identical requested baselines (`GRAD_6`, `CONE_30`, `B_10`,
`FREQ_30`, `MU_030`, `ZETA_050`, `CPAR_X1`, and `KWOB_X1`) are aliases rather than
duplicate solves.

Workflow:

1. Generate `SCREENING ONLY` STEP files with `build_screening_geometry.py`.
2. Export global-0.10 mm C3D4 meshes with Abaqus/CAE `mesh_screening_steps.py`.
3. Run `prepare_screen.py` and sequential `run_screen.ps1` jobs.
4. Run `analyze_cases.py`, `render_screen.py`, `plot_metrics.py`, and `write_report.py`.
5. Complete the blank `human_review_template.csv` while reviewing the gallery/GIFs.

Raw ODB/SIM/deck/Socket telemetry and private dense NPZ files are intentionally ignored.
The committed pose/proximity tables are lightweight 10/20 us samples reconstructed from
dense RP history. Contact events separated by at most 5 us are merged; short events
followed by at least 20 us without contact are labelled `IMPACT`, sustained events with
mean axial speed at least 50 mm/s are labelled `SLIDING`, and the remaining sustained
events are labelled `STUCK`. These labels are coarse review aids, not production contact
identification.

