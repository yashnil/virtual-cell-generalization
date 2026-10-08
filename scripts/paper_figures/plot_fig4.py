# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Figure 4: source choice strongly influences perturbation transfer (N5; N6 limitation in panel D)."""

from __future__ import annotations

import numpy as np
import style as st
from _common import load, save


def panel_forest(ax, comp):
    c = comp.sort_values("C_S", ascending=True).reset_index(drop=True)
    for i, r in c.iterrows():
        gw = r.source.startswith("K562_GWPS")
        col = st.EMPH if gw else st.ZEROSHOT
        mfc = "white" if r.source == "K562_GWPS_depth" else col
        st.herrorbar(
            ax, [r.C_S], [i], [r.lo], [r.hi], col, marker="o", mfc=mfc, ms=st.MS + 0.6, lw=1.0
        )
        ax.text(
            1.07,
            i,
            f"{r.median_source_reliability:.2f}",
            fontsize=st.FS_SMALL,
            va="center",
            ha="center",
            transform=ax.get_yaxis_transform(),
            color=st.INK if gw else st.INK2,
        )
        ax.text(
            r.lo - 0.02,
            i,
            f"{r.C_S:.2f}",
            fontsize=st.FS_SMALL,
            va="center",
            ha="right",
            color=st.INK2,
        )
    ax.text(
        1.07,
        len(c) - 0.45,
        "median\nsource\nreliability",
        fontsize=st.FS_SMALL,
        va="bottom",
        ha="center",
        transform=ax.get_yaxis_transform(),
        color=st.INK2,
        linespacing=1.0,
    )
    ax.set_yticks(range(len(c)))
    labels = []
    for r in c.itertuples():
        labels.append(r.label.replace(", depth-matched", "\n(depth-matched)"))
    ax.set_yticklabels(labels)
    for tl, r in zip(ax.get_yticklabels(), c.itertuples(), strict=True):
        if r.source.startswith("K562_GWPS"):
            tl.set_fontweight("bold")
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, len(c) - 0.4)
    st.ref_line(ax, 0, orient="v")
    ax.set_xlabel("Noise-corrected latent cosine with K562-essential\ntarget, $C_S$ (95 % CI)")
    ax.set_title("Source → K562 transfer similarity", loc="left")


SHORT = {
    "Primary: pooled latent cosine, Δ = C_GWPS − C_RPE1": "Primary: Δ$C_S$",
    "GWPS subsampled to RPE1 cell counts (Adj-4)": "GWPS at RPE1 depth (Adj-4)",
    "Raw median per-perturbation r difference": "Raw median difference",
    "Reliability + magnitude adjusted, perturbation FE (Adj-2)": "Reliability + magnitude\n+ perturbation FE (Adj-2)",
    "Adj-2 with depth-matched GWPS": "Adj-2, depth-matched",
    "Reliability-matched perturbations, n = 108 (Adj-3)": "Reliability-matched\npairs, n = 108 (Adj-3)",
}


def panel_robustness(ax, rob):
    blocks = [
        ("cosine", "Pooled latent cosine difference"),
        ("pearson", "Per-perturbation Pearson r difference"),
    ]
    y = 0
    yt, yl = [], []
    heads = []
    for key, head in blocks[::-1]:
        x = rob[rob.estimand == key].iloc[::-1]
        for r in x.itertuples():
            st.herrorbar(
                ax, [r.estimate], [y], [r.lo], [r.hi], st.EMPH, marker="o", ms=st.MS, lw=1.0
            )
            ax.text(
                r.hi + 0.012,
                y,
                f"+{r.estimate:.3f}",
                fontsize=st.FS_SMALL,
                va="center",
                color=st.INK2,
            )
            yt.append(y)
            yl.append(SHORT[r.variant])
            y += 1
        heads.append((y - 0.3, head))
        y += 0.9
    for yy, h in heads:
        ax.text(
            -0.02,
            yy,
            h,
            fontsize=st.FS_SMALL,
            fontweight="bold",
            va="bottom",
            ha="right",
            color=st.INK,
            transform=ax.get_yaxis_transform(),
        )
    ax.set_yticks(yt)
    ax.set_yticklabels(yl, fontsize=st.FS_SMALL, linespacing=1.0)
    ax.tick_params(axis="y", length=0)
    st.ref_line(ax, 0, orient="v", color=st.INK2)
    ax.set_xlim(-0.02, 0.56)
    ax.set_ylim(-0.6, y - 0.4)
    ax.set_xlabel("GWPS advantage over RPE1 (95 % CI)")
    ax.set_title("Advantage under preregistered controls", loc="left")


def panel_ecdf(ax, pp):
    for col, ls, lab, lw in (
        ("diff_full", "-", "full-depth GWPS", st.LW + 0.3),
        ("diff_depth_matched", (0, (3, 1.5)), "depth-matched GWPS", st.LW),
    ):
        v = np.sort(pp[col].to_numpy())
        f = np.arange(1, len(v) + 1) / len(v)
        fav = (v > 0).mean()
        med = np.median(v)
        ax.step(
            v,
            f,
            where="post",
            color=st.EMPH,
            ls=ls,
            lw=lw,
            label=f"{lab}\n{fav * 100:.1f} % favour GWPS\nmedian +{med:.3f}",
        )
        ax.plot(
            [med],
            [0.5],
            "o",
            color=st.EMPH,
            mfc="white" if "depth" in col else st.EMPH,
            ms=st.MS,
            zorder=4,
        )
    st.ref_line(ax, 0, orient="v", color=st.INK2)
    st.ref_line(ax, 0.5, ls=(0, (1, 2)))
    ax.set_xlim(-0.62, 0.85)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel("$r_{GWPS} − r_{RPE1}$ per perturbation (template-removed)")
    ax.set_ylabel("Cumulative fraction of perturbations")
    ax.text(
        0.83,
        0.04,
        "n = 1,054 perturbations\n○ ● medians",
        fontsize=st.FS_SMALL,
        va="bottom",
        ha="right",
        color=st.INK2,
    )
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, 0.9),
        handlelength=1.8,
        labelspacing=0.7,
        frameon=True,
        facecolor="white",
        edgecolor="none",
        framealpha=1,
    )
    ax.set_title("Advantage is broad across perturbations", loc="left")


def panel_design(fig, rect, comp):
    """D: 2 × 2 identification grid (same cell line × same study). Geometry in millimetres."""
    ax, W, H = st.schematic_axes(fig, rect)
    C = comp.set_index("source").C_S
    fs, fss = st.FS, st.FS_SMALL
    cw, chh, g = 30.0, 14.5, 1.6
    x0, y0 = 25.0, 13.0
    cols = [x0, x0 + cw + g]
    rows = [y0 + chh + g, y0]  # top = same study
    ax.text(0.0, H - 1.0, "Identification limit", fontsize=fs, va="top")
    ax.text(
        x0 + cw + g / 2,
        H - 6.0,
        "Same cell line as target (K562)?",
        fontsize=fss,
        fontweight="bold",
        ha="center",
        va="top",
    )
    for cx, lab in zip(cols, ("yes", "no"), strict=True):
        ax.text(
            cx + cw / 2,
            rows[0] + chh + 1.2,
            lab,
            fontsize=fss,
            ha="center",
            va="bottom",
            color=st.INK2,
        )
    for ry, lab in zip(
        rows, ("Same study\nas target\n(Replogle 2022)", "Different\nstudy"), strict=True
    ):
        ax.text(x0 - 1.8, ry + chh / 2, lab, fontsize=fss, ha="right", va="center", linespacing=1.1)
    ax.text(
        x0 - 13.0,
        rows[0] + chh + 1.2,
        "Same study?",
        fontsize=fss,
        fontweight="bold",
        ha="center",
        va="bottom",
    )

    def cell(
        ci, ri, title, lines, *, ec=st.INK2, lw=st.STROKE, fc="white", hatch=None, ls="-", tc=st.INK
    ):
        x, y = cols[ci], rows[ri]
        st.box(ax, x, y, cw, chh, fc=fc, ec=ec, lw=lw, hatch=hatch, ls=ls)
        if hatch:
            st.box(ax, x + 1.6, y + 1.6, cw - 3.2, chh - 3.2, fc=fc, ec="none", lw=0, zorder=1.5)
        ax.text(
            x + cw / 2,
            y + chh - 2.0,
            title,
            fontsize=fss,
            fontweight="bold",
            ha="center",
            va="top",
            color=tc,
            zorder=4,
        )
        ax.text(
            x + cw / 2,
            y + chh / 2 - 1.6,
            lines,
            fontsize=fss,
            ha="center",
            va="center",
            color=tc,
            linespacing=1.2,
            zorder=4,
        )

    cell(0, 0, "K562 GWPS", f"$C_S$ = {C['K562_GWPS']:.2f}\n(different screen)", ec=st.EMPH, lw=1.2)
    cell(1, 0, "RPE1", f"$C_S$ = {C['RPE1']:.2f}\n(same library)")
    cell(
        1,
        1,
        "Four other lines",
        f"Jurkat {C['Jurkat']:.2f}, HepG2 {C['HepG2']:.2f}\n(same lab)  ·  HCT116 {C['HCT116']:.2f},\nHEK293T {C['HEK293T']:.2f} (other lab)",
    )
    cell(
        0,
        1,
        "VIPerturb-seq K562",
        "reliability gate failed\n(0.087 < 0.10; N6)\nnot interpretable",
        fc="#F2F2F2",
        hatch="////",
        ls=(0, (3, 1.5)),
        tc=st.INK2,
    )
    ax.text(
        W / 2 + 6,
        6.0,
        "The only usable same-cell source shares the target's study and lab:\ncell identity and study/lab effects are not separable.",
        fontsize=fss,
        ha="center",
        va="center",
        style="italic",
        linespacing=1.2,
    )
    return ax


def main() -> None:
    comp = load("fig4_compatibility")
    rob = load("fig4_robustness")
    pp = load("fig4_per_perturbation")
    fig = st.new_figure(st.DOUBLE, 4.9)
    ax_a = fig.add_axes([0.14, 0.6, 0.3, 0.35])
    ax_b = fig.add_axes([0.73, 0.6, 0.235, 0.35])
    ax_c = fig.add_axes([0.085, 0.09, 0.36, 0.38])
    panel_forest(ax_a, comp)
    panel_robustness(ax_b, rob)
    panel_ecdf(ax_c, pp)
    ax_d = panel_design(fig, [0.505, 0.03, 0.485, 0.455], comp)
    st.panel_label(ax_a, "A", x=-0.115)
    st.panel_label(ax_b, "B", x=-0.215)
    st.panel_label(ax_c, "C", x=-0.06)
    st.panel_label(ax_d, "D", x=0.0, y=0.995)
    save(
        fig,
        "fig4_source_compatibility",
        sources=["fig4_compatibility", "fig4_robustness", "fig4_per_perturbation"],
        plot_script="scripts/paper_figures/plot_fig4.py",
        upstream=[
            "outputs/n5 (N5, reports/n5_results.md)",
            "reports/n6_results.md (panel D facts only)",
        ],
    )


if __name__ == "__main__":
    main()
