from __future__ import annotations

import argparse
from pathlib import Path


def node_id(i, j, k, nx, ny):
    return 1 + i + (nx + 1) * (j + (ny + 1) * k)


def elem_id(i, j, k, nx, ny):
    return 1 + i + nx * (j + ny * k)


def write_case(out_dir: Path, name: str, nonreflecting: bool, nx: int = 44, ny: int = 20, nz: int = 20):
    out_dir.mkdir(parents=True, exist_ok=True)
    dx = 12.0 / nx
    dy = dz = 0.075
    y0 = z0 = -0.75
    lines = []
    a = lines.append
    def data_chunks(values, width=16):
        vals = list(values)
        for start in range(0, len(vals), width):
            a(", ".join(str(v) for v in vals[start:start+width]))
    a("*Heading")
    a("** F100 G2P20 CEL axial boundary reflection harness")
    a("** H0: default free Eulerian axial ends; H1: documented nonreflecting axial ends")
    a("*Preprint, echo=NO, model=NO, history=NO, contact=NO")
    a("*Part, name=FLUID_EULERIAN")
    a("*Node")
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                x = -6.0 + i * dx
                y = y0 + j * dy
                z = z0 + k * dz
                a(f"{node_id(i,j,k,nx,ny)}, {x:.12g}, {y:.12g}, {z:.12g}")
    a("*Element, type=EC3D8R")
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                n000 = node_id(i, j, k, nx, ny)
                n100 = node_id(i + 1, j, k, nx, ny)
                n110 = node_id(i + 1, j + 1, k, nx, ny)
                n010 = node_id(i, j + 1, k, nx, ny)
                n001 = node_id(i, j, k + 1, nx, ny)
                n101 = node_id(i + 1, j, k + 1, nx, ny)
                n111 = node_id(i + 1, j + 1, k + 1, nx, ny)
                n011 = node_id(i, j + 1, k + 1, nx, ny)
                a(f"{elem_id(i,j,k,nx,ny)}, {n000}, {n100}, {n110}, {n010}, {n001}, {n101}, {n111}, {n011}")
    a("*Elset, elset=FLUID_ALL, generate")
    a(f"1, {nx*ny*nz}, 1")
    for label, lo, hi in [("PROBE_CENTER", 21, 23), ("PROBE_LOW_QUARTER", 10, 12), ("PROBE_HIGH_QUARTER", 32, 34), ("PROBE_LOW_END_NEAR", 2, 4), ("PROBE_HIGH_END_NEAR", 41, 43)]:
        a(f"*Elset, elset={label}")
        for k in range(nz):
            for j in range(ny):
                data_chunks(elem_id(i,j,k,nx,ny) for i in range(lo-1, hi))
        a(f"*Nset, nset={label}_NODES")
        for k in range(nz + 1):
            for j in range(ny + 1):
                data_chunks(node_id(i,j,k,nx,ny) for i in range(lo-1, hi+1))
    a("*Surface, type=ELEMENT, name=EUL_LOW")
    for k in range(nz):
        for j in range(ny):
            a(f"{elem_id(0,j,k,nx,ny)}, S1")
    a("*Surface, type=ELEMENT, name=EUL_HIGH")
    for k in range(nz):
        for j in range(ny):
            a(f"{elem_id(nx-1,j,k,nx,ny)}, S2")
    a("*Eulerian Section, elset=FLUID_ALL")
    a("MAT_FLUID_CEL, WATER")
    a("*End Part")

    a("*Part, name=CONFINING_WALL")
    a("*Node")
    wall_nodes = []
    label = 1
    surfaces = []
    # Four longitudinal planes, with one R3D4 per axial/cross-section cell.
    for side in range(4):
        grid = []
        for q in range((nz if side < 2 else ny) + 1):
            row = []
            for i in range(nx + 1):
                if side == 0:
                    xyz = (-6.0 + i*dx, y0, z0 + q*dz)
                elif side == 1:
                    xyz = (-6.0 + i*dx, y0 + ny*dy, z0 + q*dz)
                elif side == 2:
                    xyz = (-6.0 + i*dx, y0 + q*dy, z0)
                else:
                    xyz = (-6.0 + i*dx, y0 + q*dy, z0 + nz*dz)
                a(f"{label}, {xyz[0]:.12g}, {xyz[1]:.12g}, {xyz[2]:.12g}")
                row.append(label)
                label += 1
            grid.append(row)
        surfaces.append([])
        for q in range(len(grid)-1):
            for i in range(nx):
                n1, n2 = grid[q][i], grid[q][i+1]
                n3, n4 = grid[q+1][i+1], grid[q+1][i]
                el = len(wall_nodes) + 1
                wall_nodes.append((el,n1,n2,n3,n4))
                surfaces[-1].append(el)
    a("*Element, type=R3D4")
    for el,n1,n2,n3,n4 in wall_nodes:
        a(f"{el}, {n1}, {n2}, {n3}, {n4}")
    a("*Elset, elset=WALL_ALL, generate")
    a(f"1, {len(wall_nodes)}, 1")
    for idx, els in enumerate(surfaces, 1):
        a(f"*Elset, elset=WALL_SIDE_{idx}")
        data_chunks(els)
    a("*Surface, type=ELEMENT, name=WALL_SURF")
    for el, *_ in wall_nodes:
        a(f"{el}, S1")
    a("*End Part")

    a("*Assembly, name=Assembly")
    a("*Instance, name=Fluid_EULERIAN-1, part=FLUID_EULERIAN")
    a("*End Instance")
    a("*Instance, name=Confining_Wall-1, part=CONFINING_WALL")
    a("*End Instance")
    a("*Node")
    a("900000, 0., 0., 0.")
    a("*Nset, nset=RP_WALL")
    a("900000")
    a("*Elset, elset=WALL_ALL, instance=Confining_Wall-1, generate")
    a(f"1, {len(wall_nodes)}, 1")
    a("*Surface, type=EULERIAN MATERIAL, name=FLUID_SURF")
    a("Fluid_EULERIAN-1_WATER")
    a("*Rigid Body, ref node=RP_WALL, elset=WALL_ALL")
    a("*End Assembly")

    a("*Initial Conditions, type=VOLUME FRACTION")
    a("Fluid_EULERIAN-1.FLUID_ALL, Fluid_EULERIAN-1_WATER, 1.0")
    a("*Nset, nset=SRC_PLUS, instance=Fluid_EULERIAN-1")
    for k in range(nz+1):
        for j in range(ny+1):
            for i in range(21, 23):
                a(str(node_id(i,j,k,nx,ny)))
    a("*Nset, nset=SRC_MINUS, instance=Fluid_EULERIAN-1")
    for k in range(nz+1):
        for j in range(ny+1):
            for i in range(23, 25):
                a(str(node_id(i,j,k,nx,ny)))
    a("*Initial Conditions, type=VELOCITY")
    a("SRC_PLUS, 1, 10.0")
    a("SRC_MINUS, 1, -10.0")

    a("*Material, name=MAT_FLUID_CEL")
    a("*Density")
    a("1.0007e-9")
    a("*Eos, type=USUP")
    a("100000.0, 0.0, 0.0")
    a("*Viscosity")
    a("7.1e-10")
    a("*Material, name=MAT_WALL")
    a("*Density")
    a("7.8e-9")
    a("*Elastic")
    a("2.1e5, 0.3")
    a("*Surface Interaction, name=FLUID_WALL_PROP")
    a("*Surface Behavior, pressure-overclosure=HARD")
    a("*Contact")
    a("*Contact Inclusions")
    a("Confining_Wall-1.WALL_SURF, FLUID_SURF")
    a("*Contact Property Assignment")
    a("Confining_Wall-1.WALL_SURF, FLUID_SURF, FLUID_WALL_PROP")
    a("*Boundary")
    a("RP_WALL, ENCASTRE")
    a("*Step, name=Wave, nlgeom=YES")
    a("*Dynamic, Explicit")
    a(", 0.00025")
    if nonreflecting:
        a("*Eulerian Boundary, OUTFLOW=NONREFLECTING")
        a("Fluid_EULERIAN-1.EUL_LOW")
        a("Fluid_EULERIAN-1.EUL_HIGH")
    a("*Bulk Viscosity")
    a("0.06, 1.2")
    a("*Output, field, time interval=2.5e-6, time marks=YES")
    for probe in ["PROBE_CENTER", "PROBE_LOW_QUARTER", "PROBE_HIGH_QUARTER", "PROBE_LOW_END_NEAR", "PROBE_HIGH_END_NEAR"]:
        a(f"*Element Output, elset=Fluid_Eulerian-1.{probe}")
        a("EVF, S")
        a(f"*Node Output, nset=Fluid_Eulerian-1.{probe}_NODES")
        a("V")
    a("*End Step")
    (out_dir / f"{name}.inp").write_text("\n".join(lines) + "\n", encoding="ascii")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    write_case(out, "HARNESS_DEFAULT_FREE", False)
    write_case(out, "HARNESS_NONREFLECTING", True)


if __name__ == "__main__":
    main()
