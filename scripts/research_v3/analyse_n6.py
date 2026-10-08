"""Apply the preregistered N6 gate, contrasts and outcomes (n6_protocol.md §5).

Reproduce: ``uv run python scripts/research_v3/analyse_n6.py`` (after run_n6.py).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "outputs" / "n6"
CODE = [
    "reports/n6_protocol.md",
    "scripts/research_v3/build_n6.py",
    "scripts/research_v3/run_n6.py",
    "scripts/research_v3/analyse_n6.py",
    "src/virtual_cell/analysis/source_compatibility.py",
]
MIN_POOLED_REL = 0.10
MAX_DEPTH_SHIFT = 0.10
MAX_CTRL_CI_WIDTH = 0.15
MIN_PANEL = 300


def ci(x) -> dict:
    x = np.asarray(x)
    return {
        "median": float(np.median(x)),
        "lo": float(np.percentile(x, 2.5)),
        "hi": float(np.percentile(x, 97.5)),
    }


def classify(C: dict, boot: pd.DataFrame, gate: dict, adj2_vip_lo: float) -> dict:
    """Mechanical N6 outcome from point estimates, paired bootstrap and the gate."""
    d1 = gate["viperturb_pooled_split_half_reliability"] < MIN_POOLED_REL
    d2 = gate["viperturb_reliable_energy_boot_lo"] <= 0
    ctrl = ci(boot["K562_GWPS_vipdepth"])
    d3 = (
        abs(C["K562_GWPS_vipdepth"] - C["K562_GWPS"]) > MAX_DEPTH_SHIFT
        or (ctrl["hi"] - ctrl["lo"]) > MAX_CTRL_CI_WIDTH
    )
    d4 = gate["n_panel"] < MIN_PANEL
    gate_fail = {
        "d1_low_pooled_reliability": bool(d1),
        "d2_reliable_energy_not_positive": bool(d2),
        "d3_positive_control_unstable_at_vip_depth": bool(d3),
        "d4_panel_too_small": bool(d4),
    }
    vj = boot["VIPerturb_K562"] - boot["Jurkat"]
    gj = boot["K562_GWPS"] - boot["Jurkat"]
    f = vj / gj
    res = {
        "gate_failures": gate_fail,
        "VIP_minus_Jurkat": ci(vj),
        "GWPS_minus_Jurkat": ci(gj),
        "fraction_f": ci(f),
        "VIP_minus_GWPS": ci(boot["VIPerturb_K562"] - boot["K562_GWPS"]),
        "VIP_minus_RPE1": ci(boot["VIPerturb_K562"] - boot["RPE1"]),
        "positive_control_vipdepth": ctrl,
    }
    if any(gate_fail.values()):
        outcome = "D"
    elif ci(gj)["lo"] <= 0:
        outcome = "inconclusive (reference not separated from Jurkat)"
    else:
        lo_vj = ci(vj)["lo"]
        f_lo, f_hi = ci(f)["lo"], ci(f)["hi"]
        if f_hi <= 0.33:
            outcome = "B"
        elif lo_vj > 0 and adj2_vip_lo > 0 and f_lo >= 0.67:
            outcome = "A"
        elif lo_vj > 0 and adj2_vip_lo > 0:
            outcome = "C"
        else:
            outcome = "inconclusive (not robust to reliability adjustment or interval spans B/C)"
    res["outcome"] = outcome
    return res


def main() -> None:
    s = json.loads((OUT / "n6_summary.json").read_text())
    gate = json.loads((OUT / "source_gate.json").read_text())
    boot = pd.read_csv(OUT / "bootstrap_C.csv")
    adj = s["Adj2_vs_Jurkat"]["VIPerturb_K562"]
    res = classify(s["C"], boot, gate, adj["lo"])
    bb = pd.read_csv(OUT / "bootstrap_C_panelB.csv")
    res["panel_B_secondary"] = {
        "n": s["panel_B"]["n"],
        "C": s["panel_B"]["C"],
        "VIP_minus_Jurkat": ci(bb["VIPerturb_K562"] - bb["Jurkat"]),
        "fraction_f": ci((bb["VIPerturb_K562"] - bb["Jurkat"]) / (bb["K562_GWPS"] - bb["Jurkat"])),
    }
    res.update(
        {
            "C_S": {k: {"point": s["C"][k], **ci(boot[k])} for k in s["C"]},
            "Adj2_vs_Jurkat": s["Adj2_vs_Jurkat"],
            "Adj3_VIP_vs_Jurkat": s["Adj3_VIP_vs_Jurkat"],
            "secondary": {k: s[k] for k in ("R1_median_r", "R2_raw_cosine", "R3_dir_acc")},
            "gate": gate,
            "code_sha256": {f: hashlib.sha256((REPO / f).read_bytes()).hexdigest() for f in CODE},
        }
    )
    (OUT / "n6_decision.json").write_text(json.dumps(res, indent=2))
    print(
        json.dumps(
            {k: res[k] for k in ("outcome", "gate_failures", "VIP_minus_Jurkat", "fraction_f")},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
