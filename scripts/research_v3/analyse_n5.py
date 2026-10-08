"""Apply the preregistered N5 contrasts, ladder and outcome criteria (n5_protocol.md §5).

Reproduce: ``uv run python scripts/research_v3/analyse_n5.py`` (after run_n5.py).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "outputs" / "n5"
NON_K562 = ["RPE1", "HepG2", "Jurkat", "HCT116", "HEK293T"]
CODE = [
    "reports/n5_protocol.md",
    "scripts/research_v3/build_n5_gwps.py",
    "scripts/research_v3/run_n5.py",
    "scripts/research_v3/analyse_n5.py",
    "src/virtual_cell/analysis/source_compatibility.py",
]


def ci(x: np.ndarray) -> dict:
    return {
        "point_median": float(np.median(x)),
        "lo": float(np.percentile(x, 2.5)),
        "hi": float(np.percentile(x, 97.5)),
    }


def main() -> None:
    s = json.loads((OUT / "n5_summary.json").read_text())
    comp = json.loads((OUT / "matched_comparator.json").read_text())
    b = pd.read_csv(OUT / "bootstrap_C.csv")
    C = s["C"]
    contrasts = {
        "delta_BC": {"point": C["K562_GWPS"] - C["RPE1"], **ci(b.K562_GWPS - b.RPE1)},
        "delta_BC_depth": {
            "point": C["K562_GWPS_depth"] - C["RPE1"],
            **ci(b.K562_GWPS_depth - b.RPE1),
        },
        "GWPS_minus_matched": {
            "matched": comp["matched"],
            "point": C["K562_GWPS"] - C[comp["matched"]],
            **ci(b.K562_GWPS - b[comp["matched"]]),
        },
        "GWPS_minus_mean_nonK562": {
            "point": C["K562_GWPS"] - np.mean([C[x] for x in NON_K562]),
            **ci(b.K562_GWPS - b[NON_K562].mean(axis=1)),
        },
    }
    lad_b = [
        b.K562_GWPS,
        b.RPE1,
        b[["HepG2", "Jurkat"]].mean(axis=1),
        b[["HCT116", "HEK293T"]].mean(axis=1),
    ]
    names = ["GWPS", "RPE1", "mean(HepG2,Jurkat)", "mean(HCT116,HEK293T)"]
    steps = []
    for i in range(3):
        d = lad_b[i] - lad_b[i + 1]
        steps.append(
            {
                "step": f"{names[i]} > {names[i + 1]}",
                **ci(d),
                "holds": bool(np.percentile(d, 2.5) > 0),
            }
        )
    C_ci = {k: ci(b[k]) for k in b.columns}

    adj2 = s["Adj2_full"]["K562_GWPS"]
    adj2_d = s["Adj2_depth"]["K562_GWPS"]
    d_bc, d_bc_d = contrasts["delta_BC"], contrasts["delta_BC_depth"]
    A = d_bc["lo"] > 0 and adj2["lo"] > 0 and d_bc_d["lo"] > 0
    B = (s["R1_diff_GWPS_RPE1"]["lo"] > 0) and not A
    Cp = d_bc["hi"] < 0
    Dp = C_ci["K562_GWPS"]["hi"] < 0.70
    spread = max(C[k] for k in ["K562_GWPS"] + NON_K562) - min(
        C[k] for k in ["K562_GWPS"] + NON_K562
    )
    E = spread < 0.05 and d_bc["lo"] <= 0 <= d_bc["hi"]
    outcomes = {
        "A_same_cell_advantage_survives": bool(A),
        "B_disappears_after_control": bool(B),
        "C_prime_screen_library_dominates": bool(Cp),
        "D_prime_same_cell_mismatch_large": bool(Dp),
        "E_no_structure": bool(E),
    }
    if not (A or B or Cp or E):
        outcomes["inconclusive"] = True
    decision = {
        "C_S": {k: {"point": C[k], **C_ci[k]} for k in C},
        "contrasts": contrasts,
        "ladder": steps,
        "ladder_holds": all(x["holds"] for x in steps),
        "Adj2_GWPS_minus_RPE1": {"full": adj2, "depth": adj2_d},
        "Adj3": s["Adj3"],
        "R1_diff_GWPS_RPE1": s["R1_diff_GWPS_RPE1"],
        "spread_C": spread,
        "outcomes": outcomes,
        "kill_direction": bool(E or B or Cp),
        "code_sha256": {f: hashlib.sha256((REPO / f).read_bytes()).hexdigest() for f in CODE},
    }
    r4 = OUT / "r4_calibration.csv"
    if r4.exists():
        df = pd.read_csv(r4)
        decision["R4"] = (
            df.groupby(["source", "k"])[["M3_own", "M1"]]
            .median()
            .round(4)
            .reset_index()
            .to_dict("records")
        )
    (OUT / "n5_decision.json").write_text(json.dumps(decision, indent=2))
    print(
        json.dumps(
            {k: decision[k] for k in ("outcomes", "kill_direction", "ladder_holds")}, indent=2
        )
    )


if __name__ == "__main__":
    main()
