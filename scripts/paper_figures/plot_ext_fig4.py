# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Extended Data Fig. 4: the N6 reliability gate and its positive control (Outcome D).

VIPerturb-seq compatibility values are deliberately not shown: under the frozen protocol they are
non-interpretable after the gate failure (reports/n6_results.md §2).
"""

from __future__ import annotations

import numpy as np
import style as st
from _common import load, save

SHORT = {
    "VIPerturb_K562": "VIPerturb-seq",
    "K562_GWPS_vipdepth": "GWPS, VIPerturb depth",
    "K562_GWPS": "GWPS, full depth",
}
LAB = {
    "VIPerturb_K562": "VIPerturb-seq K562",
    "K562_GWPS_vipdepth": "K562 GWPS at VIPerturb depth",
    "K562_GWPS": "K562 GWPS, full depth",
}


def panel_pooled(ax, p):
    p = p.iloc[::-1].reset_index(drop=True)
    for i, r in p.iterrows():
        frozen_stat = r.status.startswith("frozen")
        col = st.FAILED if frozen_stat else st.CONTROL
        mfc = col if frozen_stat else "white"
        ax.plot(
            r.pooled_split_half_reliability,
            i,
            "o" if frozen_stat else "s",
            color=st.INK2 if frozen_stat else st.CONTROL,
            mfc=mfc,
            ms=st.MS + 1,
        )
        ax.text(
            r.pooled_split_half_reliability + 0.012,
            i + 0.22,
            f"{r.pooled_split_half_reliability:.4f}"
            if frozen_stat
            else f"{r.pooled_split_half_reliability:.3f}"
            + ("" if frozen_stat else "  (exploratory)"),
            fontsize=st.FS_SMALL,
            va="bottom",
            color=st.INK,
        )
    ax.axvline(0.10, color=st.INK, lw=0.9)
    ax.text(0.103, -0.5, "preregistered gate 0.10", fontsize=st.FS_SMALL, color=st.INK, va="bottom")
    ax.set_yticks(range(len(p)))
    ax.set_yticklabels(
        [f"{s}\n(median {int(c)} cells)" for s, c in zip(p.source, p.median_cells, strict=True)],
        linespacing=1.0,
    )
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 0.4)
    ax.set_ylim(-0.6, len(p) - 0.3)
    ax.set_xlabel("Pooled split-half reliability")
    ax.set_title("Gate d1 failed: 0.0866 < 0.10", loc="left")


def panel_ecdf(ax, t):
    styles = {
        "VIPerturb_K562": ("-", st.INK2, 1.2),
        "K562_GWPS_vipdepth": ((0, (3, 1.5)), st.CONTROL, 1.0),
        "K562_GWPS": (":", st.CONTROL, 1.0),
    }
    for s_, (ls, col, lw) in styles.items():
        v = np.sort(t[t.source == s_].rel.to_numpy())
        f = np.arange(1, len(v) + 1) / len(v)
        ax.step(
            v,
            f,
            where="post",
            ls=ls,
            color=col,
            lw=lw,
            label=f"{SHORT[s_]}: {np.mean(v >= 0.10) * 100:.0f} %",
        )
    ax.axvline(0.10, color=st.INK, lw=0.7, ls=(0, (2, 2)))
    st.ref_line(ax, 0, orient="v")
    ax.set_xlim(-0.15, 0.8)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel("Per-perturbation split-half reliability (n = 637)")
    ax.set_ylabel("Cumulative fraction")
    ax.legend(
        loc="lower right",
        fontsize=st.FS_SMALL,
        handlelength=2.2,
        title="≥ 0.10 per perturbation",
        title_fontsize=st.FS_SMALL,
    )
    ax.set_title("Per-perturbation reliability", loc="left")


def panel_control(ax, pc):
    full = pc.iloc[0]
    ax.axvspan(full.C_S - 0.10, full.C_S + 0.10, color="#EFEFEF", lw=0, zorder=0)
    ax.text(
        full.C_S - 0.095,
        1.45,
        "d3 tolerance ± 0.10",
        fontsize=st.FS_SMALL,
        color=st.INK2,
        va="bottom",
    )
    for i, r in pc.iloc[::-1].reset_index(drop=True).iterrows():
        st.herrorbar(
            ax,
            [r.C_S],
            [i],
            [r.lo],
            [r.hi],
            st.CONTROL,
            marker="s",
            mfc="white" if "VIPerturb" in r.source else st.CONTROL,
            ms=st.MS + 0.6,
            lw=1.0,
        )
        ax.text(
            r.hi + 0.01,
            i,
            f"{r.C_S:.3f} [{r.lo:.3f}, {r.hi:.3f}]",
            fontsize=st.FS_SMALL,
            va="center",
        )
    ax.set_yticks([0, 1])
    ax.set_yticklabels(list(pc.source.iloc[::-1]))
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0.6, 1.0)
    ax.set_ylim(-0.6, 1.9)
    ax.set_xlabel("$C_S$ vs K562-essential target (N6 panel; 95 % CI)")
    ax.set_title("Positive control passed (d3)", loc="left")


def panel_gate(ax, g):
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.text(
        0.0,
        9.6,
        "Preregistered gate (any failure → Outcome D)",
        fontsize=st.FS_SMALL,
        fontweight="bold",
        va="top",
    )
    for i, r in enumerate(g.itertuples()):
        y = 7.9 - i * 1.6
        failed = bool(r.failed)
        ax.plot(
            0.35,
            y,
            "o",
            ms=6,
            color=st.INK if failed else st.INK2,
            mfc=st.INK if failed else "white",
            mew=0.8,
        )
        ax.text(0.85, y, f"{r.item}  {r.condition_to_pass}", fontsize=st.FS_SMALL, va="center")
        ax.text(
            9.9,
            y,
            f"{r.observed}  {'FAIL' if failed else 'pass'}",
            fontsize=st.FS_SMALL,
            va="center",
            ha="right",
            color=st.INK if failed else st.INK2,
            fontweight="bold" if failed else "normal",
        )
    ax.text(
        0.0,
        0.6,
        "VIPerturb compatibility values are omitted: under the frozen\nprotocol they cannot be read for or against same-cell transfer.",
        fontsize=st.FS_SMALL,
        color=st.INK2,
        style="italic",
        va="center",
        linespacing=1.15,
    )


def main() -> None:
    fig = st.new_figure(st.DOUBLE, 4.2)
    ax_a = fig.add_axes([0.2, 0.6, 0.27, 0.33])
    ax_b = fig.add_axes([0.6, 0.6, 0.37, 0.33])
    ax_c = fig.add_axes([0.2, 0.1, 0.25, 0.3])
    ax_d = fig.add_axes([0.55, 0.04, 0.43, 0.42])
    panel_pooled(ax_a, load("ext4_pooled_reliability"))
    panel_ecdf(ax_b, load("ext4_per_perturbation_rel"))
    pc = load("ext4_positive_control")
    panel_control(ax_c, pc)
    panel_gate(ax_d, load("ext4_gate"))
    st.panel_label(ax_a, "A", x=-0.17)
    st.panel_label(ax_b, "B", x=-0.07)
    st.panel_label(ax_c, "C", x=-0.17)
    st.panel_label(ax_d, "D", x=0.0, y=0.98)
    save(
        fig,
        "ext_fig4_n6_reliability",
        sources=[
            "ext4_pooled_reliability",
            "ext4_per_perturbation_rel",
            "ext4_positive_control",
            "ext4_gate",
        ],
        plot_script="scripts/paper_figures/plot_ext_fig4.py",
        upstream=["outputs/n6", "reports/n6_results.md §4 (exploratory values)"],
    )


if __name__ == "__main__":
    main()
