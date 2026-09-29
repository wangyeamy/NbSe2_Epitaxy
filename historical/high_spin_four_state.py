#!/usr/bin/env python3
"""Analyze only certified k3 high-spin four-state endpoints."""
from __future__ import annotations
import argparse, csv, hashlib, json, math, re
from pathlib import Path

ROOT = Path("./REDACTED_PATH")
MECH = Path("./REDACTED_PATH")
OLD_C = Path("./REDACTED_PATH")
OLD_H = Path("./REDACTED_PATH")
RY_EV = 13.605693122994
BOHR_A = 0.529177210903
FORCE_CONV = RY_EV / BOHR_A
CELL = ((14.2770005232, 0.0), (-7.1385002616, 12.3642393781))

STATES = {
    "U_GU": {
        "marker": OLD_C / "results/k3_interface.converged", "nat": 183,
        "energy": -12271.36385020, "moment": 10.33, "branch": "high",
    },
    "U_GC": {
        "marker": Path("./REDACTED_PATH"),
        "nat": 183, "energy": -12271.36382170, "moment": 10.34, "branch": "high",
        "charge": "14bae5f3c1fdfc27b26d747ca5ddd47cc49eda08cc40885cefba550b7fca5b72",
    },
    "C_GU": {
        "marker": Path("./REDACTED_PATH"),
        "nat": 245, "energy": -13507.58049110, "moment": 10.41, "branch": "high",
        "charge": "7fc426f541b94da74fe0178dedcfed73bd8257e42c36627da6841602117064f3",
    },
    "C_GC": {
        "marker": Path("./REDACTED_PATH"),
        "nat": 245, "energy": -13507.58072249, "moment": 10.41, "branch": "high",
        "charge": "6d1201b7b5d329a731b71e4187f87c1d508e3dca0aaa8262c3c4dfe8cd6ad4d9",
    },
    "H": {
        "marker": OLD_H / "results/k3_hbn.converged", "nat": 62, "branch": "reference",
    },
}

def file_sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""): digest.update(block)
    return digest.hexdigest()

def marker_record(path):
    if not path.is_file(): raise FileNotFoundError(path)
    text = path.read_text()
    try:
        value = json.loads(text)
        if not isinstance(value, dict): raise ValueError("marker JSON is not an object")
        return value
    except json.JSONDecodeError:
        return {key.strip(): value.strip() for line in text.splitlines() if "=" in line for key, value in [line.split("=", 1)]}

def parse_qe(path, nat):
    if not path.is_file(): raise FileNotFoundError(path)
    text = path.read_text(errors="replace")
    if "convergence has been achieved" not in text or "JOB DONE." not in text:
        raise ValueError(f"uncertified output: {path}")
    if "convergence NOT achieved" in text or "Error in routine" in text:
        raise ValueError(f"failure marker in output: {path}")
    def last(pattern):
        values = re.findall(pattern, text, re.M)
        if not values: raise ValueError(f"missing {pattern}: {path}")
        return float(values[-1].replace("D", "E"))
    energy = last(r"^!\s+total energy\s+=\s+([-+0-9.EeDd]+)\s+Ry")
    moment = last(r"total magnetization\s+=\s+([-+0-9.EeDd]+)\s+Bohr mag/cell")
    absolute = last(r"absolute magnetization\s+=\s+([-+0-9.EeDd]+)\s+Bohr mag/cell")
    accuracy = last(r"estimated scf accuracy\s+<\s+([-+0-9.EeDd]+)\s+Ry")
    lines = text.splitlines(); blocks = []
    for index, line in enumerate(lines):
        if "Forces acting on atoms (cartesian axes, Ry/au):" not in line: continue
        block = []
        for row in lines[index+1:]:
            match = re.search(r"atom\s+\d+\s+type\s+\d+\s+force =\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)", row)
            if match: block.append(tuple(float(x.replace("D", "E")) for x in match.groups()))
            elif block: break
        if len(block) == nat: blocks.append(block)
    if not blocks: raise ValueError(f"no complete {nat}-atom force block: {path}")
    return {"output": str(path), "output_sha256": file_sha(path), "energy_ry": energy, "moment": moment, "absolute_moment": absolute, "accuracy_ry": accuracy, "forces": blocks[-1]}

def load_states():
    result = {}
    for name, spec in STATES.items():
        record = marker_record(spec["marker"])
        output = Path(record.get("output", ""))
        parsed = parse_qe(output, spec["nat"])
        if "energy" in spec and abs(parsed["energy_ry"] - spec["energy"]) > 1e-7:
            raise ValueError(f"{name} energy lock failed")
        if "moment" in spec and abs(parsed["moment"] - spec["moment"]) > 0.05:
            raise ValueError(f"{name} magnetic branch lock failed")
        if "charge" in spec and record.get("charge_sha256") != spec["charge"]:
            raise ValueError(f"{name} certified charge lock failed")
        if parsed["accuracy_ry"] >= 1e-8:
            raise ValueError(f"{name} accuracy gate failed")
        parsed.update({"marker": str(spec["marker"]), "marker_sha256": file_sha(spec["marker"]), "nat": spec["nat"], "branch": spec["branch"]})
        result[name] = parsed
    return result

def read_geometry(path, expected, expected_sha):
    if file_sha(path) != expected_sha: raise ValueError(f"geometry hash mismatch: {path}")
    with path.open(newline="") as stream:
        rows = [(row["species"], float(row["x"]), float(row["y"]), float(row["z"])) for row in csv.DictReader(stream)]
    if len(rows) != expected: raise ValueError(f"geometry atom count mismatch: {path}")
    return rows

def subtract(a, b): return [tuple(x-y for x, y in zip(va, vb)) for va, vb in zip(a, b)]
def norm(v): return math.sqrt(sum(x*x for x in v))
def metrics(vectors):
    values = [norm(row) for row in vectors]
    return {"rms_ev_a": math.sqrt(sum(x*x for x in values)/len(values))*FORCE_CONV, "max_ev_a": max(values)*FORCE_CONV, "sum_ev_a": [sum(row[j] for row in vectors)*FORCE_CONV for j in range(3)]}

def mode_vectors(atoms):
    layer = atoms[135:183]
    nbz = sum(z for symbol, x, y, z in layer if symbol == "Nb") / 16
    raw = {"shear_x": [], "shear_y": [], "breathing": [], "buckling": [], "rigid_z": []}
    for symbol, x, y, z in layer:
        top = symbol == "Se" and z > nbz; bottom = symbol == "Se" and z < nbz
        raw["shear_x"].append((1 if top else -1 if bottom else 0, 0, 0))
        raw["shear_y"].append((0, 1 if top else -1 if bottom else 0, 0))
        raw["breathing"].append((0, 0, 1 if top else -1 if bottom else 0))
        raw["buckling"].append((0, 0, 1 if symbol == "Nb" else -0.5))
        raw["rigid_z"].append((0, 0, 1))
    result = {}
    for name, rows in raw.items():
        scale = math.sqrt(sum(sum(x*x for x in row) for row in rows))
        result[name] = [tuple(x/scale for x in row) for row in rows]
    return result

def projections(vectors, atoms):
    layer = vectors[135:183]; modes = mode_vectors(atoms); result = {}; internal = []
    for name, mode in modes.items():
        value = sum(sum(force[j]*basis[j] for j in range(3)) for force, basis in zip(layer, mode))*FORCE_CONV
        result[name + "_ev_a_unit_mode"] = value
        if name != "rigid_z": internal.append(value/FORCE_CONV)
    total = sum(sum(x*x for x in row) for row in layer)
    result["internal_mode_span_fraction"] = sum(x*x for x in internal)/total if total else 0.0
    result["shear_magnitude_ev_a_unit_mode"] = math.hypot(result["shear_x_ev_a_unit_mode"], result["shear_y_ev_a_unit_mode"])
    return result

def unwrap_delta(a, b):
    dx, dy = a[0]-b[0], a[1]-b[1]
    det = CELL[0][0]*CELL[1][1]-CELL[0][1]*CELL[1][0]
    u = (dx*CELL[1][1]-dy*CELL[1][0])/det; v = (CELL[0][0]*dy-CELL[0][1]*dx)/det
    u -= round(u); v -= round(v)
    return (u*CELL[0][0]+v*CELL[1][0], u*CELL[0][1]+v*CELL[1][1], a[2]-b[2])

def run(check_only=False):
    gu = read_geometry(MECH/"source/GU.csv", 183, "3089572f80d1c94ca97da20733576c0978b52182bbef5a969be439355bec3c1e")
    gc_hc = read_geometry(MECH/"source/GC_HC.csv", 245, "7c01b6c6329310b35d025f9c011cd89b66c4c8207a8260eee8a2c372b03534c6")
    gc = gc_hc[:183]
    states = load_states()
    if check_only:
        print("HIGH_SPIN_FOUR_STATE_REFERENCES_VERIFIED")
        return
    cap_gu = subtract(states["C_GU"]["forces"][:183], states["U_GU"]["forces"])
    cap_gc = subtract(states["C_GC"]["forces"][:183], states["U_GC"]["forces"])
    nonadd = subtract(cap_gc, cap_gu)
    hbn_gu = subtract(states["C_GU"]["forces"][183:], states["H"]["forces"])
    hbn_gc = subtract(states["C_GC"]["forces"][183:], states["H"]["forces"])
    delta_u = (states["U_GC"]["energy_ry"]-states["U_GU"]["energy_ry"])*RY_EV
    delta_c = (states["C_GC"]["energy_ry"]-states["C_GU"]["energy_ry"])*RY_EV
    ecap_gu = (states["C_GU"]["energy_ry"]-states["U_GU"]["energy_ry"]-states["H"]["energy_ry"])*RY_EV
    ecap_gc = (states["C_GC"]["energy_ry"]-states["U_GC"]["energy_ry"]-states["H"]["energy_ry"])*RY_EV
    displacement = [unwrap_delta(a[1:], b[1:]) for a, b in zip(gc, gu)]
    summary = {
        "scope": "certified k3 fixed-geometry high-spin four-state analysis; no gap scan or reaction barrier",
        "energy": {"delta_U_GC_minus_GU_ev": delta_u, "delta_C_GC_minus_GU_ev": delta_c, "capping_interaction_GU_ev": ecap_gu, "capping_interaction_GC_ev": ecap_gc, "delta_delta_capping_GC_minus_GU_ev": ecap_gc-ecap_gu},
        "forces": {
            "capping_on_substrate_GU": {"metrics": metrics(cap_gu), "modes": projections(cap_gu, gu)},
            "capping_on_substrate_GC": {"metrics": metrics(cap_gc), "modes": projections(cap_gc, gc)},
            "four_state_nonadditivity": {"metrics": metrics(nonadd), "modes_on_GC_basis": projections(nonadd, gc)},
            "hbn_interaction_GU": {"metrics": metrics(hbn_gu)}, "hbn_interaction_GC": {"metrics": metrics(hbn_gc)},
        },
        "GU_to_GC_displacement_modes": {key.replace("_ev_a_unit_mode", "_angstrom_unit_mode"): value/FORCE_CONV for key, value in projections(displacement, gc).items() if key.endswith("_ev_a_unit_mode")},
        "states": {name: {key: value for key, value in state.items() if key != "forces"} for name, state in states.items()},
    }
    out = ROOT/"results"; out.mkdir(exist_ok=True)
    (out/"high_spin_four_state_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n")
    with (out/"substrate_force_differences.csv").open("w", newline="") as stream:
        fields = ["atom_index", "species", "GU_x_A", "GU_y_A", "GU_z_A", "GC_x_A", "GC_y_A", "GC_z_A", "cap_GU_fx_eV_A", "cap_GU_fy_eV_A", "cap_GU_fz_eV_A", "cap_GC_fx_eV_A", "cap_GC_fy_eV_A", "cap_GC_fz_eV_A", "nonadd_fx_eV_A", "nonadd_fy_eV_A", "nonadd_fz_eV_A"]
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for index, (g1, g2, f1, f2, fd) in enumerate(zip(gu, gc, cap_gu, cap_gc, nonadd), 1):
            writer.writerow(dict(zip(fields, [index, g1[0], *g1[1:], *g2[1:], *(x*FORCE_CONV for x in f1), *(x*FORCE_CONV for x in f2), *(x*FORCE_CONV for x in fd)])))
    marker = {"completed": True, "summary": str(out/"high_spin_four_state_summary.json"), "summary_sha256": file_sha(out/"high_spin_four_state_summary.json"), "force_table": str(out/"substrate_force_differences.csv"), "force_table_sha256": file_sha(out/"substrate_force_differences.csv"), "delta_delta_capping_meV_cell": (ecap_gc-ecap_gu)*1000}
    (out/"HIGH_SPIN_FOUR_STATE_COMPLETE.json").write_text(json.dumps(marker, indent=2, sort_keys=True)+"\n")
    print("HIGH_SPIN_FOUR_STATE_ANALYSIS_COMPLETE")
    print(json.dumps(marker, indent=2, sort_keys=True))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--check-only", action="store_true")
    run(parser.parse_args().check_only)
