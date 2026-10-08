"""Apply the preregistered N1 decision rule and draw N1/N4 figures.

Reads only outputs of run_n1_calibration_budget.py / run_n4_agreement_null.py.
Rules are transcribed from reports/n1_n4_protocol.md §2.7 and §3.5.

Reproduce: ``uv run python scripts/analyse_n1_n4.py``
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
N1 = REPO / "outputs" / "n1_n4" / "n1"
N4 = REPO / "outputs" / "n1_n4" / "n4"
FIG = REPO / "outputs" / "n1_n4" / "figures"
SRC = REPO / "data" / "figure_sources" / "n1_n4"
K_REF = 885
ORDER = ["K562", "RPE1", "HepG2", "Jurkat"]
METRICS = ["M0", "M1", "M2", "M3", "M4", "M1r", "M3r"]
CODE = [
    "src/virtual_cell/analysis/calibration_budget.py",
    "src/virtual_cell/analysis/agreement_null.py",
    "scripts/build_n1_n4_splits.py",
    "scripts/run_n1_calibration_budget.py",
    "scripts/run_n4_agreement_null.py",
    "scripts/analyse_n1_n4.py",
    "reports/n1_n4_protocol.md",
]
COLORS = {
    "E0": "#9aa0a6",
    "E0s": "#5f6368",
    "E1": "#e8a33d",
    "E2": "#c5672b",
    "E3": "#3b7dd8",
    "E4": "#7b3fbf",
}


def summarise(draws: pd.DataFrame) -> pd.DataFrame:
    g = draws.groupby(["context", "estimator", "k"])
    rows = []
    for (c, e, k), d in g:
        row = {"context": c, "estimator": e, "k": k, "n_draws": len(d)}
        for m in METRICS:
            row[f"{m}_median"] = d[m].median()
            row[f"{m}_lo"] = d[m].quantile(0.025)
            row[f"{m}_hi"] = d[m].quantile(0.975)
        rows.append(row)
    return pd.DataFrame(rows)


def decide(summary: pd.DataFrame, boot: pd.DataFrame) -> dict:
    out = {}
    b3 = boot[boot.metric == "M3"].set_index(["context", "estimator", "k"])
    for est in ("E3", "E4"):
        per = {}
        for c in ORDER:
            s = summary[(summary.context == c) & (summary.estimator == est)].set_index("k")
            k20, kref = s.loc[20], s.loc[K_REF]
            c1_draw = bool(k20["M3_lo"] > 0)
            c1_boot = bool(b3.loc[(c, est, 20), "boot_lo"] > 0)
            c2 = bool(kref["M3_median"] > 0)
            c3 = bool(c2 and k20["M3_median"] >= 0.25 * kref["M3_median"])
            ref_boot_pos = bool(b3.loc[(c, est, K_REF), "boot_lo"] > 0)
            onset = None
            for k in sorted(s.index):
                if k == 0:
                    continue
                if s.loc[k, "M3_lo"] > 0 and b3.loc[(c, est, k), "boot_lo"] > 0:
                    onset = int(k)
                    break
            per[c] = {
                "M3_k20_median": float(k20["M3_median"]),
                "M3_k20_draw_lo": float(k20["M3_lo"]),
                "M3_k20_boot_lo": float(b3.loc[(c, est, 20), "boot_lo"]),
                "M3_kref_median": float(kref["M3_median"]),
                "M3_kref_boot_lo": float(b3.loc[(c, est, K_REF), "boot_lo"]),
                "C1": c1_draw and c1_boot,
                "C2": c2,
                "C3": c3,
                "pass": c1_draw and c1_boot and c2 and c3,
                "failA_condition": c2 and ref_boot_pos,
                "gamma_onset_k": onset,
            }
        out[est] = {
            "per_context": per,
            "n_pass": sum(v["pass"] for v in per.values()),
            "n_failA": sum(v["failA_condition"] for v in per.values()),
        }
    if any(out[e]["n_pass"] >= 3 for e in ("E3", "E4")):
        verdict = "PASS"
    elif any(out[e]["n_failA"] >= 3 for e in ("E3", "E4")):
        verdict = "FAIL-A"
    else:
        verdict = "FAIL-B"
    out["N1"] = verdict
    return out


def curves(summary: pd.DataFrame, metric: str, ests, fname: str, title: str) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8), sharey=False)
    for ax, c in zip(axes, ORDER, strict=True):
        for e in ests:
            s = summary[(summary.context == c) & (summary.estimator == e)].sort_values("k")
            if s.empty:
                continue
            x = s.k.replace(0, 0.5).to_numpy()
            ax.plot(x, s[f"{metric}_median"], "-o", ms=3, color=COLORS[e], label=e)
            ax.fill_between(x, s[f"{metric}_lo"], s[f"{metric}_hi"], color=COLORS[e], alpha=0.15)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xscale("log")
        ax.set_title(c)
        ax.set_xlabel("k anchors (0 plotted at 0.5; rightmost = k_ref 885)")
    axes[0].set_ylabel(metric)
    axes[-1].legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=150)
    plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    SRC.mkdir(parents=True, exist_ok=True)
    draws = pd.read_csv(N1 / "n1_draws.csv")
    boot = pd.read_csv(N1 / "n1_boot.csv")
    summary = summarise(draws)
    summary.to_csv(N1 / "n1_summary.csv", index=False)
    summary.to_csv(SRC / "n1_summary.csv", index=False)
    decision = decide(summary, boot)

    # secondary: gamma_ANOVA gain beyond scale (paired by draw)
    piv = draws[draws.k > 0].pivot_table(
        index=["context", "k", "repeat", "draw"], columns="estimator", values="M4"
    )
    gains = []
    for e in ("E3", "E4"):
        d = (piv[e] - piv["E2"]).groupby(level=["context", "k"])
        for (c, k), v in d:
            gains.append(
                {
                    "context": c,
                    "k": k,
                    "estimator": e,
                    "M4_gain_over_E2_median": v.median(),
                    "lo": v.quantile(0.025),
                    "hi": v.quantile(0.975),
                }
            )
    pd.DataFrame(gains).to_csv(N1 / "n1_m4_gain_over_e2.csv", index=False)

    # F3a anchor-permutation null
    f3a = pd.read_csv(N1 / "n1_f3a.csv")
    f3a_s = f3a.groupby(["context", "estimator", "k"]).M3.median().rename("perm_M3_median")
    unperm = summary.set_index(["context", "estimator", "k"])["M3_median"]
    flags = []
    for (c, e, k), v in f3a_s.items():
        u = unperm.loc[(c, e, k)]
        flags.append(
            {
                "context": c,
                "estimator": e,
                "k": k,
                "permuted_M3_median": v,
                "unpermuted_M3_median": u,
                "flag": bool(u > 0 and v >= 0.5 * u),
            }
        )
    pd.DataFrame(flags).to_csv(N1 / "n1_f3a_summary.csv", index=False)
    decision["F3a_flags"] = [f for f in flags if f["flag"]]

    # F3b shared-control diagnostic (paired by draw)
    f3b = pd.read_csv(N1 / "n1_f3b.csv")
    f3b_s = f3b.groupby(["context", "estimator", "k", "variant"])[["M0", "M3"]].median().unstack()
    f3b_s.to_csv(N1 / "n1_f3b_summary.csv")

    decision["code_sha256"] = {f: hashlib.sha256((REPO / f).read_bytes()).hexdigest() for f in CODE}
    (N1 / "n1_decision.json").write_text(json.dumps(decision, indent=2, default=str))

    ests = ["E0", "E0s", "E1", "E2", "E3", "E4"]
    curves(
        summary,
        "M0",
        ests,
        "n1_M0_full_raw.png",
        "M0: full response incl. template (energy explained)",
    )
    curves(
        summary,
        "M1",
        ests,
        "n1_M1_full_template_removed.png",
        "M1: full response, template removed",
    )
    curves(summary, "M2", ests, "n1_M2_conserved.png", "M2: conserved (along source consensus)")
    curves(
        summary,
        "M3",
        ["E2", "E3", "E4"],
        "n1_M3_gamma_perp.png",
        "M3: γ⊥ (orthogonal to source consensus) — primary endpoint",
    )
    curves(summary, "M4", ests, "n1_M4_gamma_anova.png", "M4: γ_ANOVA residual removed")
    curves(summary, "M3r", ["E3", "E4"], "n1_M3r_gamma_perp_corr.png", "M3r: median r(P⊥, y⊥)")
    print(
        json.dumps({k: v for k, v in decision.items() if k != "code_sha256"}, indent=2, default=str)
    )

    if (N4 / "n4_t3_null.csv").exists():
        null = pd.read_csv(N4 / "n4_t3_null.csv")
        t3 = pd.read_csv(N4 / "n4_t3.csv").set_index("context")
        fig, axes = plt.subplots(1, 4, figsize=(14, 3.2))
        for ax, c in zip(axes, ORDER, strict=True):
            ax.hist(null[c], bins=30, color="#bbb")
            ax.axvline(t3.loc[c, "observed_centred"], color="#c5672b", lw=2)
            ax.set_title(f"{c}: {t3.loc[c, 'class']}")
            ax.set_xlabel("Spearman(agreement, −D)")
        fig.tight_layout()
        fig.savefig(FIG / "n4_t3_exchangeable_null.png", dpi=150)
        plt.close(fig)


if __name__ == "__main__":
    main()
