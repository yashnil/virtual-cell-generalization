"""Figure 11 — coverage is not transferability (C1 vs C1 + KOLF2.1J, C4).

Left: Arc targets supported by >= k direct sources. Right: held-out change from adding
KOLF as one more equal-weight donor. Drawn only from
``data/figure_sources/fig11_coverage_vs_transferability.csv``.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from virtual_cell.visualization import common, style


def main() -> None:
    style.apply()
    t = common.load_source("fig11_coverage_vs_transferability")
    cov = t[t.panel == "coverage"]
    tr = t[t.panel == "transfer"]
    fig, (ax0, ax1) = plt.subplots(
        1, 2, figsize=(style.FULL_WIDTH, 3.0), gridspec_kw={"width_ratios": [1, 1.5]}
    )
    y = np.arange(len(cov))
    ax0.barh(y + 0.2, cov.c1, 0.38, color=style.GREY, label="C1 (3 GREEN atlases)")
    ax0.barh(y - 0.2, cov.c1_kolf, 0.38, color=style.BLUE, label="C1 + KOLF2.1J")
    for i, (a, b) in enumerate(zip(cov.c1, cov.c1_kolf, strict=True)):
        ax0.text(b + 3, i - 0.2, f"{b:.0f}", va="center", fontsize=7)
        ax0.text(a + 3, i + 0.2, f"{a:.0f}", va="center", fontsize=7, color=style.GREY)
    ax0.set_yticks(y, cov.quantity, fontsize=7)
    ax0.invert_yaxis()
    ax0.set_xlim(0, 330)
    ax0.set_xlabel("Arc targets (of 300)")
    ax0.set_title("Coverage gained", fontsize=9)
    ax0.legend(fontsize=6, loc="lower right")

    labels = [f"{f}: {q}" for f, q in zip(tr.fold, tr.quantity, strict=True)]
    colors = [style.GREEN if d > 0 else style.VERMILION for d in tr.delta]
    yy = np.arange(len(tr))
    ax1.barh(yy, tr.delta, color=colors, edgecolor="black", lw=0.4)
    for i, d in enumerate(tr.delta):
        ax1.text(d, i, f" {d:+.3f}", va="center", ha="left" if d >= 0 else "right", fontsize=6)
    style.reference_line(ax1, 0, axis="x")
    ax1.set_yticks(yy, labels, fontsize=7)
    ax1.invert_yaxis()
    ax1.set_xlabel("change from adding KOLF (held-out atlas)")
    ax1.set_title("Transfer: helps only the pluripotent fold", fontsize=9)
    lo, hi = tr.delta.min(), tr.delta.max()
    pad = 0.35 * (hi - lo)
    ax1.set_xlim(lo - pad, hi + pad)
    fig.suptitle(
        "Coverage is not transferability: KOLF fills 11 of 13 unsupported Arc targets "
        "but worsens transfer outside a pluripotent context",
        fontsize=8,
    )
    fig.tight_layout()
    common.save_figure(fig, "fig11_coverage_vs_transferability")
    print("  wrote reports/figures/fig11_coverage_vs_transferability.{png,svg}")


if __name__ == "__main__":
    main()
