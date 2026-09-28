from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ZIP_PATH = ROOT / "CEL_OVERNIGHT_12H_RESULTS.zip"
MANIFEST_PATH = ROOT / "CEL_OVERNIGHT_12H_RESULTS_MANIFEST.json"
EXCLUSIONS_PATH = ROOT / "CEL_OVERNIGHT_12H_RESULTS_EXCLUSIONS.txt"

ROOT_FILES = {
    "CEL_OVERNIGHT_12H_ROOTCAUSE.md",
    "OVERNIGHT_CEL_STATE.json",
    "CEL_FLUID_CONTACT_FRICTION_comparison.gif",
    "F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX20_robot_motion.gif",
    "F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX_FLUIDFRIC0_20_robot_motion.gif",
    "F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX_RF_FRIC0_20_robot_motion.gif",
    "F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX_WF_FRIC0_20_robot_motion.gif",
    "analyze_initialization_closure_results.py",
    "analyze_fluidfric0_closure.py",
    "analyze_rf_fric0_closure.py",
    "analyze_wf_fric0_closure.py",
    "audit_initfix_final.py",
    "extract_initfix_motion.py",
    "prepare_initialization_fix.py",
    "render_friction_comparison_gif.py",
}

CASE_DIRS = {
    "F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX20",
    "F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX_FLUIDFRIC0_20",
    "F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX_RF_FRIC0_20",
    "F100_G2P20_ZEROPRESSURE_FACEFIX_INITFIX_WF_FRIC0_20",
}
CASE_EXTENSIONS = {".inp", ".f", ".for", ".f90", ".dat", ".sta", ".log", ".json", ".csv"}
CLOSURE_EXTENSIONS = {".json", ".csv", ".md", ".txt", ".log"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def select_files() -> list[Path]:
    selected: set[Path] = {ROOT / name for name in ROOT_FILES if (ROOT / name).is_file()}
    closure = ROOT / "initialization_closure"
    for path in closure.rglob("*"):
        if path.is_file() and path.suffix.lower() in CLOSURE_EXTENSIONS:
            selected.add(path)
    cases = ROOT / "case"
    for case_name in CASE_DIRS:
        case_dir = cases / case_name
        for path in case_dir.rglob("*"):
            if path.is_file() and path.suffix.lower() in CASE_EXTENSIONS:
                selected.add(path)
    return sorted(path for path in selected if path.is_file())


def main() -> None:
    files = select_files()
    missing = sorted(name for name in ROOT_FILES if not (ROOT / name).is_file())
    entries = []
    for path in files:
        rel = path.relative_to(ROOT).as_posix()
        entries.append({"path": rel, "bytes": path.stat().st_size, "sha256": sha256(path)})

    excluded_lines = [
        "Excluded by design from CEL_OVERNIGHT_12H_RESULTS.zip:",
        "- all *.odb, *.sim, *.simdir, *.res, *.pac, *.prt, *.sel, *.stt, *.abq, *.mdl, *.msg, *.lck, *.env files",
        "- all restart files and private *.npz histories",
        "- unrelated cases and unrelated working-tree files",
        "",
        f"Selected files: {len(entries)}",
        f"Missing expected root artifacts: {', '.join(missing) if missing else 'none'}",
    ]
    EXCLUSIONS_PATH.write_text("\n".join(excluded_lines) + "\n", encoding="utf-8")

    manifest = {
        "package": ZIP_PATH.name,
        "generated_from": str(ROOT),
        "root_cause_report": "CEL_OVERNIGHT_12H_ROOTCAUSE.md",
        "selected_file_count": len(entries),
        "excluded_artifact_policy": "ODB/SIM/RES/restart/NPZ and unrelated files excluded",
        "files": entries,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in files:
            archive.write(path, path.relative_to(ROOT).as_posix())
        archive.write(MANIFEST_PATH, MANIFEST_PATH.name)
        archive.write(EXCLUSIONS_PATH, EXCLUSIONS_PATH.name)
    print(json.dumps({"zip": str(ZIP_PATH), "bytes": ZIP_PATH.stat().st_size, "files": len(entries)}, indent=2))


if __name__ == "__main__":
    main()
