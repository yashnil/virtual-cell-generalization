"""Figure 9 — does G1 control transport deliver the response it is given while
keeping single-cell structure? (frozen summary statistics only)"""

from __future__ import annotations

import matplotlib.pyplot as plt

from virtual_cell.visualization import common, style

GROUPS = [
    ("real_control_reference", "real controls (reference)", style.BLACK, "o"),
    ("real_perturbed", "real perturbed cells", style.GREEN, "D"),
    ("G0_control_resample", "G0 control resampling", style.GREY, "s"),
    ("G1_transport", "G1 transport (s = 0.5)", style.BLUE, "o"),
    ("G1_transport_unsmoothed", "G1 unsmoothed (s = 0)", style.SKY, "v"),
    ("G2_count_model", "G2 count model", style.ORANGE, "^"),
]


def main() -> None:
    style.apply()
    t = common.load_source("fig9_transport_fidelity")
    st = t[t.panel == "structure"]
    fi = t[t.panel == "fidelity"]
    fig, axes = plt.subplots(2, 2, figsize=(style.FULL_WIDTH, 5.2))

    ax = axes[0, 0]
    names = ["G1_transport", "G1_transport_unsmoothed", "G2_count_model", "G0_control_resample"]
    labels = {g: lab for g, lab, _, _ in GROUPS}
    for i, g in enumerate(names):
        rec = fi[fi.group == g].set_index("stat").value
        color = {n: c for n, _, c, _ in GROUPS}[g]
        ax.barh(i - 0.18, rec["realised_norm"], height=0.34, color=color, edgecolor="black", lw=0.6)
        ax.barh(
            i + 0.18,
            rec["own_intended_norm"],
            height=0.34,
            color="white",
            edgecolor=color,
            hatch="//",
            lw=0.8,
        )
        txt = (
            f"slope {rec['slope_vs_own']:.3f}, r {rec['pearson_vs_own']:.3f}"
            if rec["own_intended_norm"] > 0
            else "asked for no effect"
        )
        ax.text(
            max(rec["realised_norm"], rec["own_intended_norm"]) + 0.8,
            i,
            txt,
            va="center",
            fontsize=6.8,
        )
    g0 = float(fi[(fi.group == "G0_control_resample") & (fi.stat == "realised_norm")].value.iloc[0])
    style.reference_line(ax, g0, axis="x", color=style.GREY)
    ax.set_yticks(range(len(names)), [labels[g] for g in names], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 55)
    ax.set_xlabel(
        "pseudobulk response norm (50 perturbations)\n"
        "solid = realised, hatched = intended,\n"
        "dashed = 400-cell noise floor (G0)"
    )
    ax.set_title("Mean-effect delivery", fontsize=9)
    style.panel_label(ax, "A")

    def dot_range(ax, quantity, title, letter, xlabel):
        for i, (g, _lab, color, marker) in enumerate(GROUPS):
            rec = st[(st.group == g) & (st.quantity == quantity)].set_index("stat").value
            ax.plot([rec["q10"], rec["q90"]], [i, i], color=color, lw=2)
            ax.plot(rec["median"], i, marker=marker, color=color, mec="black", ms=6, ls="")
        ax.set_yticks(range(len(GROUPS)), [lab for _, lab, _, _ in GROUPS], fontsize=7)
        ax.invert_yaxis()
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(xlabel)
        style.panel_label(ax, letter)

    dot_range(
        axes[0, 1], "library_size", "Library size per cell", "B", "UMIs per cell (median, q10–q90)"
    )
    axes[0, 1].set_xlim(0, 25_000)
    dot_range(
        axes[1, 0],
        "genes_detected",
        "Genes detected per cell",
        "C",
        "genes detected (median, q10–q90)",
    )
    axes[1, 0].set_xlim(0, 5_000)

    ax = axes[1, 1]
    for i, (g, _lab, color, _marker) in enumerate(GROUPS):
        sp = float(st[(st.group == g) & (st.quantity == "sparsity")].value.iloc[0])
        ax.barh(
            i,
            sp,
            color=color,
            edgecolor="black",
            lw=0.6,
            height=0.6,
            hatch="" if g.startswith("real") else "//",
        )
        ax.text(sp + 0.01, i, f"{sp:.3f}", va="center", fontsize=7)
    ref = float(st[(st.group == "real_perturbed") & (st.quantity == "sparsity")].value.iloc[0])
    style.reference_line(ax, ref, axis="x", color=style.GREEN)
    ax.set_yticks(range(len(GROUPS)), [lab for _, lab, _, _ in GROUPS], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("zero fraction (sparsity)")
    ax.set_title("Sparsity", fontsize=9)
    style.panel_label(ax, "D")
    fig.text(
        0.01,
        0.005,
        "Held-out K562, 300 perturbations × 400 generated cells. Only frozen "
        "summaries are shown (quantiles, fitted slope/r); per-gene and per-cell values were "
        "not retained.",
        fontsize=6.3,
        color=style.GREY,
    )
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    common.save_figure(fig, "fig9_transport_fidelity")


if __name__ == "__main__":
    main()
