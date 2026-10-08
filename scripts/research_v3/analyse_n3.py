"""Apply the preregistered N3 rules (reports/n3_protocol.md) and draw figures.

Reads only N3 outputs plus source-only/basal quantities from the N3 data object.

Reproduce: ``uv run python scripts/research_v3/analyse_n3.py``
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
N3 = REPO / "outputs" / "n3"
DATA = N3 / "data"
FIG = N3 / "figures"
SRC = REPO / "data" / "figure_sources" / "n3"
CONTEXTS = ["K562", "RPE1", "HepG2", "Jurkat", "HCT116", "HEK293T"]
ORIG, NEW = CONTEXTS[:4], CONTEXTS[4:]
K_REF = 743
GRID = [5, 10, 20, 50, 100, 200, K_REF]
SEED = 20261006
METRICS = ["M0", "M1", "M2", "M3", "M4", "M1r", "M3r"]
COLORS = {
    "E0": "#9aa0a6",
    "E0s": "#5f6368",
    "E1": "#e8a33d",
    "E2": "#c5672b",
    "E3": "#3b7dd8",
    "E4": "#7b3fbf",
}
CODE = [
    "scripts/research_v3/build_n3_six_context.py",
    "scripts/research_v3/run_n3a.py",
    "scripts/research_v3/run_n3b.py",
    "scripts/research_v3/run_n3c.py",
    "scripts/research_v3/run_n3_n4.py",
    "scripts/research_v3/analyse_n3.py",
    "src/virtual_cell/analysis/calibration_budget.py",
    "src/virtual_cell/analysis/agreement_null.py",
    "reports/n3_protocol.md",
]


def q(x, p):
    return float(np.percentile(x, p))


# --------------------------------------------------------------------------
# N3-A
# --------------------------------------------------------------------------


def summarise(draws: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (c, e, k), d in draws.groupby(["context", "estimator", "k"]):
        row = {"context": c, "estimator": e, "k": k, "n_draws": len(d)}
        for m in METRICS:
            row[f"{m}_median"] = d[m].median()
            row[f"{m}_lo"] = d[m].quantile(0.025)
            row[f"{m}_hi"] = d[m].quantile(0.975)
        rows.append(row)
    return pd.DataFrame(rows)


def first_k(series: pd.Series, threshold: float):
    for k in sorted(series.index):
        if k >= 1 and series.loc[k] >= threshold:
            return int(k)
    return None


def n3a_decide(s: pd.DataFrame, boot: pd.DataFrame) -> dict:
    b3 = boot[boot.metric == "M3"].set_index(["context", "estimator", "k"])
    out = {}
    for c in CONTEXTS:
        per = {}
        for est in ("E3", "E4"):
            x = s[(s.context == c) & (s.estimator == est)].set_index("k")
            k20, kref = x.loc[20], x.loc[K_REF]
            c1 = bool(k20["M3_lo"] > 0 and b3.loc[(c, est, 20), "boot_lo"] > 0)
            c2 = bool(kref["M3_median"] > 0)
            c3 = bool(c2 and k20["M3_median"] >= 0.25 * kref["M3_median"])
            per[est] = {
                "C1": c1,
                "C2": c2,
                "C3": c3,
                "pass": c1 and c2 and c3,
                "failA_cond": c2 and bool(b3.loc[(c, est, K_REF), "boot_lo"] > 0),
                "M3_k20": float(k20["M3_median"]),
                "M3_kref": float(kref["M3_median"]),
                "M3_kref_boot_lo": float(b3.loc[(c, est, K_REF), "boot_lo"]),
            }
        if any(v["pass"] for v in per.values()):
            cls = "PASS"
        elif any(v["failA_cond"] for v in per.values()):
            cls = "FAIL-A-type"
        else:
            cls = "FAIL-B-type"
        # pattern components
        e2 = s[(s.context == c) & (s.estimator == "E2")].set_index("k")["M0_median"]
        e0s = float(s[(s.context == c) & (s.estimator == "E0s") & (s.k == 0)]["M0_median"].iloc[0])
        gain = e2 - e0s
        kT = first_k(gain, 0.5 * gain.loc[K_REF]) if gain.loc[K_REF] > 0 else None
        m3 = s[(s.context == c) & (s.estimator == "E4")].set_index("k")["M3_median"]
        kG = first_k(m3, 0.5 * m3.loc[K_REF]) if m3.loc[K_REF] > 0 else None
        out[c] = {
            "class": cls,
            "estimators": per,
            "k_T50": kT if gain.loc[K_REF] > 0 else "no template/scale gain",
            "k_gamma50": kG if m3.loc[K_REF] > 0 else "none",
            "P1_template_small_k": bool(kT is not None and kT <= 20),
            "P2_gamma_slower": bool(kT is not None and kG is not None and kG > kT),
            "P3_low_budget_minority": bool(m3.loc[K_REF] > 0 and m3.loc[20] / m3.loc[K_REF] < 0.25),
            "E4_M3_ratio_k20_kref": float(m3.loc[20] / m3.loc[K_REF])
            if m3.loc[K_REF] > 0
            else None,
            "E2_M0_gain_kref": float(gain.loc[K_REF]),
            "E0s_M1_k0": float(
                s[(s.context == c) & (s.estimator == "E0s") & (s.k == 0)]["M1_median"].iloc[0]
            ),
        }
    classes = {c: out[c]["class"] for c in NEW}
    out["R_A1_failA_replicates_externally"] = all(v == "FAIL-A-type" for v in classes.values())
    return out


def curves(s, metric, ests, fname, title):
    fig, axes = plt.subplots(1, 6, figsize=(22, 3.8))
    for ax, c in zip(axes, CONTEXTS, strict=True):
        for e in ests:
            x = s[(s.context == c) & (s.estimator == e)].sort_values("k")
            if x.empty:
                continue
            kk = x.k.replace(0, 0.5).to_numpy()
            ax.plot(kk, x[f"{metric}_median"], "-o", ms=3, color=COLORS[e], label=e)
            ax.fill_between(kk, x[f"{metric}_lo"], x[f"{metric}_hi"], color=COLORS[e], alpha=0.15)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xscale("log")
        ax.set_title(c + (" (X-Atlas)" if c in NEW else ""))
        ax.set_xlabel("k (rightmost = k_ref 743)")
    axes[0].set_ylabel(metric)
    axes[-1].legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=140)
    plt.close(fig)


# --------------------------------------------------------------------------
# N3-B
# --------------------------------------------------------------------------


def kstar(med: pd.DataFrame, col: int, gstar: float) -> float:
    """Smallest budget with median R >= g*, log-linear interpolation between grid points."""
    ks = [k for k in GRID if k in med.index]
    vals = [med.loc[k, col] for k in ks]
    if vals[0] >= gstar:
        return float(ks[0])
    for k0, v0, k1, v1 in zip(ks, vals, ks[1:], vals[1:], strict=False):
        if v0 < gstar <= v1:
            f = (gstar - v0) / (v1 - v0)
            return float(np.exp(np.log(k0) + f * (np.log(k1) - np.log(k0))))
    return float("inf")


def n3b(dr: pd.DataFrame, n3a_draws: pd.DataFrame) -> dict:
    e4 = dr[dr.estimator == "E4"]
    # R(k, m): per draw mean over subsets of size m, then median over draws
    per_draw = e4.groupby(["context", "k", "m", "repeat", "draw"])[
        ["M3_ref", "M3_own", "M1"]
    ].mean()
    R = (
        per_draw.groupby(["context", "k", "m"])
        .agg(
            R_gamma=("M3_ref", "median"),
            R_gamma_lo=("M3_ref", lambda x: q(x, 2.5)),
            R_gamma_hi=("M3_ref", lambda x: q(x, 97.5)),
            R_gamma_own=("M3_own", "median"),
            R_full=("M1", "median"),
            R_full_lo=("M1", lambda x: q(x, 2.5)),
            R_full_hi=("M1", lambda x: q(x, 97.5)),
        )
        .reset_index()
    )
    e0 = (
        dr[(dr.estimator == "E0s") & (dr.k == 0)]
        .groupby(["context", "m", "repeat"])[["M3_ref", "M1"]]
        .mean()
        .groupby(["context", "m"])
        .median()
        .reset_index()
    )
    e0.columns = ["context", "m", "R_gamma", "R_full"]
    e0["k"] = 0
    R = pd.concat([R, e0], ignore_index=True)

    # F-B1: m = 5 equals N3-A E4 on the same draws
    a = n3a_draws[
        (n3a_draws.estimator == "E4")
        & n3a_draws.k.isin([5, 10, 20, 50, 100, 200])
        & (n3a_draws.draw < 8)
    ].set_index(["context", "k", "repeat", "draw"])["M3"]
    b = e4[e4.m == 5].set_index(["context", "k", "repeat", "draw"])["M3_ref"]
    j = a.to_frame().join(b, how="inner")
    fb1 = float(np.max(np.abs(j.M3 - j.M3_ref))) if len(j) else float("nan")

    res = {"F_B1_max_abs_diff": fb1, "n_F_B1": int(len(j)), "per_target": {}}
    pd_ = per_draw.reset_index()
    for c in CONTEXTS:
        x = pd_[pd_.context == c]
        piv = x.pivot_table(index=["k", "repeat", "draw"], columns="m", values="M3_ref")
        d20 = piv.loc[20][5] - piv.loc[20][3]
        d20_1 = piv.loc[20][5] - piv.loc[20][1]
        r20_3 = float(piv.loc[20][3].median())
        material = bool(
            d20.median() >= 0.03
            and q(d20, 2.5) > 0
            and (r20_3 <= 0 or d20.median() >= 0.25 * r20_3)
        )
        med = R[(R.context == c) & (R.k > 0)].pivot_table(index="k", columns="m", values="R_gamma")
        gstar = float(med.loc[50, 3])

        kst = {m: kstar(med, m, gstar) for m in range(1, 6)}
        rho = kst[3] / kst[5] if np.isfinite(kst[5]) and kst[5] > 0 else float("nan")
        # combination check at k = 50: full set vs best single source, per draw
        singles = e4[(e4.context == c) & (e4.k == 50) & (e4.m == 1)].pivot_table(
            index=["repeat", "draw"], columns="subset", values="M3_ref"
        )
        full5 = e4[(e4.context == c) & (e4.k == 50) & (e4.m == 5)].set_index(["repeat", "draw"])[
            "M3_ref"
        ]
        comb = full5 - singles.max(axis=1)
        res["per_target"][c] = {
            "delta20_5v3_median": float(d20.median()),
            "delta20_5v3_lo": q(d20, 2.5),
            "delta20_5v3_hi": q(d20, 97.5),
            "delta20_5v1_median": float(d20_1.median()),
            "R20_m3": r20_3,
            "R20_m5": float(piv.loc[20][5].median()),
            "Q1_material": material,
            "g_star_R50_m3": gstar,
            "k_star": {str(m): kst[m] for m in kst},
            "rho_k3_over_k5": rho,
            "Q2_material": bool(np.isfinite(rho) and rho >= 1.5),
            "combination_k50_median": float(comb.median()),
            "combination_k50_lo": q(comb, 2.5),
            "combination_k50_hi": q(comb, 97.5),
        }
    n_q1 = sum(v["Q1_material"] for v in res["per_target"].values())
    n_q2 = sum(v["Q2_material"] for v in res["per_target"].values())
    res["n_Q1_material"], res["n_Q2_material"] = n_q1, n_q2
    res["outcome_B_criterion"] = bool(n_q1 >= 4 and n_q2 >= 4)
    return res, R


def descriptors() -> tuple[dict, np.ndarray, np.ndarray]:
    pert = np.load(DATA / "pert_part_means.npy", mmap_mode="r")
    ctrl = np.load(DATA / "ctrl_part_means.npy")
    D = np.load(DATA / "delta6.npy").astype(np.float64)
    basal = np.load(DATA / "control_means6.npy").astype(np.float64)
    rel = {}
    for c in range(6):
        vals = []
        for r in range(pert.shape[0]):
            h1 = np.asarray(pert[r, 0, c], dtype=np.float64) - ctrl[r, 0, c]
            h2 = (
                (np.asarray(pert[r, 1, c], dtype=np.float64) - ctrl[r, 1, c])
                + (np.asarray(pert[r, 2, c], dtype=np.float64) - ctrl[r, 2, c])
            ) / 2
            h1c, h2c = h1 - h1.mean(axis=1, keepdims=True), h2 - h2.mean(axis=1, keepdims=True)
            vals.append(
                np.einsum("pg,pg->p", h1c, h2c)
                / (np.linalg.norm(h1c, axis=1) * np.linalg.norm(h2c, axis=1))
            )
        rel[c] = float(np.median(np.mean(vals, axis=0)))
    Dc = D - D.mean(axis=1, keepdims=True)
    rcor = np.zeros((6, 6))
    for a in range(6):
        for b in range(6):
            x = Dc[a] - Dc[a].mean(axis=1, keepdims=True)
            y = Dc[b] - Dc[b].mean(axis=1, keepdims=True)
            rcor[a, b] = np.median(
                np.einsum("pg,pg->p", x, y)
                / (np.linalg.norm(x, axis=1) * np.linalg.norm(y, axis=1))
            )
    bsim = np.corrcoef(basal)
    return rel, rcor, bsim


def q3(dr: pd.DataFrame, rel: dict, rcor: np.ndarray, bsim: np.ndarray, k: int) -> dict:
    idx = {c: i for i, c in enumerate(CONTEXTS)}
    rho = {c: 2 * rel[c] / (1 + rel[c]) for c in range(6)}
    e4 = dr[(dr.estimator == "E4") & (dr.k == k) & dr.m.isin([2, 3, 4])]
    Rs = e4.groupby(["context", "m", "subset"]).M3_ref.median().reset_index()
    rows = []
    for _, row in Rs.iterrows():
        mem = [idx[x] for x in row.subset.split("+")]
        t = idx[row.context]
        pairs = [(a, b) for i, a in enumerate(mem) for b in mem[i + 1 :]]
        rows.append(
            {
                "context": row.context,
                "m": row.m,
                "R": row.M3_ref,
                "Rel": np.mean([rel[s] for s in mem]),
                "Div": np.mean([1 - rcor[a, b] / np.sqrt(rho[a] * rho[b]) for a, b in pairs]),
                "Partner": max(bsim[t, s] for s in mem),
            }
        )
    df = pd.DataFrame(rows)

    def zrank(g):
        out = g.copy()
        for col in ("R", "Rel", "Div", "Partner"):
            rk = stats.rankdata(g[col])
            out[col] = (rk - rk.mean()) / (rk.std() + 1e-12)
        return out

    z = df.groupby(["context", "m"], group_keys=False).apply(zrank)

    def partial(y, x, covs):
        C = np.column_stack(covs)
        ry = y - C @ np.linalg.lstsq(C, y, rcond=None)[0]
        rx = x - C @ np.linalg.lstsq(C, x, rcond=None)[0]
        return float(np.corrcoef(ry, rx)[0, 1])

    obs = {
        "Div|Rel": partial(z.R.to_numpy(), z.Div.to_numpy(), [z.Rel.to_numpy()]),
        "Rel|Div": partial(z.R.to_numpy(), z.Rel.to_numpy(), [z.Div.to_numpy()]),
        "Partner|Rel,Div": partial(
            z.R.to_numpy(), z.Partner.to_numpy(), [z.Rel.to_numpy(), z.Div.to_numpy()]
        ),
    }
    rng = np.random.default_rng([SEED, 17])
    blocks = [
        np.flatnonzero(((z.context == c) & (z.m == m)).to_numpy())
        for c, m in z[["context", "m"]].drop_duplicates().itertuples(index=False)
    ]
    null = {key: [] for key in obs}
    for _ in range(2000):
        Rp = z.R.to_numpy().copy()
        for b in blocks:
            Rp[b] = Rp[rng.permutation(b)]
        null["Div|Rel"].append(partial(Rp, z.Div.to_numpy(), [z.Rel.to_numpy()]))
        null["Rel|Div"].append(partial(Rp, z.Rel.to_numpy(), [z.Div.to_numpy()]))
        null["Partner|Rel,Div"].append(
            partial(Rp, z.Partner.to_numpy(), [z.Rel.to_numpy(), z.Div.to_numpy()])
        )
    p = {key: float((np.sum(np.abs(null[key]) >= abs(obs[key])) + 1) / 2001) for key in obs}
    df.to_csv(N3 / "n3b" / f"n3b_q3_subsets_k{k}.csv", index=False)
    return {
        "k": k,
        "partial": obs,
        "p_two_sided": p,
        "n_points": int(len(df)),
        "complementarity": bool(obs["Div|Rel"] > 0 and p["Div|Rel"] < 0.05),
        "averaging_reliability": bool(obs["Rel|Div"] > 0 and p["Rel|Div"] < 0.05),
    }


# --------------------------------------------------------------------------
# N3-C
# --------------------------------------------------------------------------


def n3c(est: pd.DataFrame, stab: pd.DataFrame, data: pd.DataFrame, n3a: dict) -> dict:
    out = {}
    for c in CONTEXTS:
        x = est[est.context == c]
        per_k = {}
        for k, g in x.groupby("k"):
            full = g.M3_full.median()
            f = {
                r: (g[f"M3_rank{r}"].median() / full if full > 0.01 else None) for r in range(1, 51)
            }
            reff = next((r for r in range(1, 51) if f[r] is not None and f[r] >= 0.9), None)
            per_k[int(k)] = {
                "M3_full_median": float(full),
                "f": {r: f[r] for r in (1, 2, 3, 5, 10)},
                "cum_gain": {r: float(g[f"M3_rank{r}"].median()) for r in (1, 2, 3, 5, 10)},
                "r_eff90": reff if full > 0.01 else "no gain",
                "participation_ratio": float(g.pr.median()),
            }
        st = stab[stab.context == c].set_index("k")
        s20 = st.loc[20]
        stable20 = bool(s20.v1_abscos_median > s20.null_v1_abscos_p975)
        f1_20 = per_k[20]["f"][1]
        f1_ref = per_k[K_REF]["f"][1]
        r20, rref = per_k[20]["r_eff90"], per_k[K_REF]["r_eff90"]
        kT = n3a[c]["k_T50"]
        h = {
            "h1": bool(isinstance(kT, int) and kT <= 20),
            "h2": bool(f1_20 is not None and f1_20 >= 0.40),
            "h3": bool(f1_20 is not None and f1_ref is not None and f1_ref < f1_20),
            "h4": bool(
                isinstance(r20, int) and (rref is None or (isinstance(rref, int) and rref > r20))
            ),
        }
        out[c] = {
            "per_k": per_k,
            "stability": st.reset_index().to_dict("records"),
            "axis_stable_vs_null_k20": stable20,
            "hierarchy": h,
            "hierarchy_reproduces": all(h.values()),
            "data_side_o_r": data[data.context == c].set_index("rank").o_r.to_dict(),
        }
    out["outcome_D_criterion"] = all(
        out[c]["hierarchy"]["h2"] and out[c]["axis_stable_vs_null_k20"] for c in NEW
    )
    return out


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    SRC.mkdir(parents=True, exist_ok=True)
    draws = pd.read_csv(N3 / "n3a" / "n3a_draws.csv")
    boot = pd.read_csv(N3 / "n3a" / "n3a_boot.csv")
    s = summarise(draws)
    s.to_csv(N3 / "n3a" / "n3a_summary.csv", index=False)
    s.to_csv(SRC / "n3a_summary.csv", index=False)
    A = n3a_decide(s, boot)

    f3a = pd.read_csv(N3 / "n3a" / "n3a_f3a.csv").groupby(["context", "estimator", "k"]).M3.median()
    un = s.set_index(["context", "estimator", "k"]).M3_median
    A["F3a_flags"] = [
        {
            "context": c,
            "estimator": e,
            "k": int(k),
            "perm": float(v),
            "unperm": float(un.loc[(c, e, k)]),
        }
        for (c, e, k), v in f3a.items()
        if un.loc[(c, e, k)] > 0 and v >= 0.5 * un.loc[(c, e, k)]
    ]
    f3b = (
        pd.read_csv(N3 / "n3a" / "n3a_f3b.csv")
        .groupby(["context", "estimator", "k", "variant"])[["M0", "M3"]]
        .median()
        .unstack()
    )
    f3b.to_csv(N3 / "n3a" / "n3a_f3b_summary.csv")

    dr = pd.read_csv(N3 / "n3b" / "n3b_draws.csv")
    B, R = n3b(dr, draws)
    R.to_csv(N3 / "n3b" / "n3b_R_gamma.csv", index=False)
    R.to_csv(SRC / "n3b_R_gamma.csv", index=False)
    rel, rcor, bsim = descriptors()
    B["Q3_k50"] = q3(dr, rel, rcor, bsim, 50)
    B["Q3_k20"] = q3(dr, rel, rcor, bsim, 20)

    C = n3c(
        pd.read_csv(N3 / "n3c" / "n3c_est.csv"),
        pd.read_csv(N3 / "n3c" / "n3c_stab.csv"),
        pd.read_csv(N3 / "n3c" / "n3c_data.csv"),
        A,
    )
    N4 = json.loads((N3 / "n4" / "n4_decision.json").read_text())

    # assay diagnostics
    D = np.load(DATA / "delta6.npy").astype(np.float64)
    genes = (REPO / "data/splits/six_context_n3/shared_genes.txt").read_text().split()
    perts = (REPO / "data/splits/six_context_n3/shared_perturbations.txt").read_text().split()
    gi = {g: i for i, g in enumerate(genes)}
    on = [(i, gi[p]) for i, p in enumerate(perts) if p in gi]
    diag = {
        CONTEXTS[c]: {
            "median_norm_delta": float(np.median(np.linalg.norm(D[c], axis=1))),
            "median_reliability": rel[c],
            "median_on_target_delta": float(np.median([D[c, i, g] for i, g in on])),
        }
        for c in range(6)
    }
    diag["basal_similarity"] = (
        pd.DataFrame(bsim, index=CONTEXTS, columns=CONTEXTS).round(3).to_dict()
    )
    diag["response_correlation"] = (
        pd.DataFrame(rcor, index=CONTEXTS, columns=CONTEXTS).round(3).to_dict()
    )

    # outcomes
    orig_failA = sum(A[c]["class"] == "FAIL-A-type" for c in ORIG)
    m3ref = {c: A[c]["estimators"]["E4"]["M3_kref"] for c in CONTEXTS}
    c1 = any(A[c]["class"] != "FAIL-A-type" for c in NEW) and orig_failA >= 3
    c2 = any(
        m3ref[c] < 0.5 * min(m3ref[o] for o in ORIG) or m3ref[c] > 2 * max(m3ref[o] for o in ORIG)
        for c in NEW
    )
    c3 = any(A[c]["E0s_M1_k0"] <= 0 for c in NEW) and all(A[o]["E0s_M1_k0"] > 0 for o in ORIG)
    outcomes = {
        "A": bool(A["R_A1_failA_replicates_externally"] and not B["outcome_B_criterion"]),
        "B": B["outcome_B_criterion"],
        "C": bool(c1 or c2 or c3),
        "C_parts": {"c1": bool(c1), "c2": bool(c2), "c3": bool(c3)},
        "D": C["outcome_D_criterion"],
    }
    any_pattern = any(
        A[c]["P1_template_small_k"] or A[c]["P2_gamma_slower"] or A[c]["P3_low_budget_minority"]
        for c in NEW
    )
    if any(A[c]["class"] == "PASS" for c in NEW) or outcomes["B"]:
        direction = "retain"
    elif outcomes["A"]:
        direction = "revise"
    elif (all(A[c]["class"] == "FAIL-B-type" for c in NEW) and not outcomes["B"]) or (
        (c1 or c2) and not any_pattern
    ):
        direction = "kill"
    else:
        direction = "unresolved"
    decision = {
        "N3A": A,
        "N3B": B,
        "N3C": C,
        "N4_six_fold": N4,
        "assay_diagnostics": diag,
        "outcomes": outcomes,
        "direction": direction,
        "code_sha256": {f: hashlib.sha256((REPO / f).read_bytes()).hexdigest() for f in CODE},
    }
    (N3 / "n3_decision.json").write_text(json.dumps(decision, indent=2, default=str))

    curves(
        s,
        "M0",
        ["E0s", "E1", "E2", "E3", "E4"],
        "n3a_M0.png",
        "N3-A M0: full response incl. template",
    )
    curves(s, "M1", ["E0s", "E2", "E3", "E4"], "n3a_M1.png", "N3-A M1: template removed")
    curves(s, "M3", ["E2", "E3", "E4"], "n3a_M3.png", "N3-A M3: γ⊥ (primary)")
    fig, axes = plt.subplots(1, 6, figsize=(22, 3.8))
    cm = plt.get_cmap("viridis")
    for ax, c in zip(axes, CONTEXTS, strict=True):
        for m in range(1, 6):
            x = R[(R.context == c) & (R.m == m) & (R.k > 0)].sort_values("k")
            ax.plot(x.k, x.R_gamma, "-o", ms=3, color=cm((m - 1) / 4), label=f"m={m}")
            ax.fill_between(x.k, x.R_gamma_lo, x.R_gamma_hi, color=cm((m - 1) / 4), alpha=0.1)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xscale("log")
        ax.set_title(c)
        ax.set_xlabel("k")
    axes[0].set_ylabel("R_γ(k, m)  (E4, M3_ref)")
    axes[-1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "n3b_R_gamma_k_m.png", dpi=140)
    plt.close(fig)
    print(
        json.dumps(
            {
                "N3A": {c: A[c]["class"] for c in CONTEXTS},
                "R_A1": A["R_A1_failA_replicates_externally"],
                "N3B": {
                    k: B[k]
                    for k in (
                        "n_Q1_material",
                        "n_Q2_material",
                        "outcome_B_criterion",
                        "F_B1_max_abs_diff",
                    )
                },
                "Q3_k50": B["Q3_k50"],
                "N3C_D": C["outcome_D_criterion"],
                "N4": N4["N4"],
                "outcomes": outcomes,
                "direction": direction,
            },
            indent=2,
            default=str,
        )
    )


if __name__ == "__main__":
    main()
