# ruff: noqa: E501 (figure labels and provenance notes are kept as single strings)
"""Automated layout QC for the paper figures: font sizes, text clipping and text-text overlaps.

Rebuilds each figure through its plot module (``save`` is intercepted, nothing is written), then inspects every
visible text artist. Writes ``reports/paper_figures/figure_qc_auto.json``.
"""

from __future__ import annotations

import importlib
import itertools
import json

import matplotlib

matplotlib.use("Agg")
import _common  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

MODULES = {
    "fig1_transfer_components": "plot_fig1",
    "fig2_sample_complexity": "plot_fig2",
    "fig3_source_count": "plot_fig3",
    "fig4_source_compatibility": "plot_fig4",
    "ext_fig1_robustness": "plot_ext_fig1",
    "ext_fig2_n1_controls": "plot_ext_fig2",
    "ext_fig3_agreement_correction": "plot_ext_fig3",
    "ext_fig4_n6_reliability": "plot_ext_fig4",
    "ext_fig5_n3_reference_frame": "plot_ext_fig5",
}


def inspect(fig) -> dict:
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    owner = {}
    for ax in fig.axes:
        for t in ax.get_xticklabels(which="both"):
            owner[id(t)] = (ax, "x")
        for t in ax.get_yticklabels(which="both"):
            owner[id(t)] = (ax, "y")

    def shown(t):
        if not (t.get_visible() and t.get_text().strip()):
            return False
        if id(t) in owner:
            ax, which = owner[id(t)]
            if not ax.axison or not (ax.xaxis if which == "x" else ax.yaxis).get_visible():
                return False
            lo, hi = sorted(ax.get_xlim() if which == "x" else ax.get_ylim())
            v = t.get_position()[0] if which == "x" else t.get_position()[1]
            return lo - 1e-9 <= v <= hi + 1e-9
        return True

    texts = [t for t in fig.findobj(matplotlib.text.Text) if shown(t)]
    sizes = [t.get_fontsize() for t in texts]
    # text box only (for annotations, exclude the leader line)
    boxes = [(t, matplotlib.text.Text.get_window_extent(t, renderer=r)) for t in texts]
    clipped = [
        t.get_text()[:40]
        for t, b in boxes
        if b.x0 < -1 or b.y0 < -1 or b.x1 > W + 1 or b.y1 > H + 1
    ]
    overlaps = []
    for (t1, b1), (t2, b2) in itertools.combinations(boxes, 2):
        if b1.overlaps(b2):
            ix = min(b1.x1, b2.x1) - max(b1.x0, b2.x0)
            iy = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)
            if ix > 2 and iy > 2:
                overlaps.append(
                    [
                        t1.get_text()[:30],
                        t2.get_text()[:30],
                        round(ix * iy / min(b1.width * b1.height, b2.width * b2.height), 2),
                    ]
                )
    # text sitting on data: data vertices / markers (display coords) inside any text box
    pts = []
    for ax in fig.axes:
        if not ax.axison:
            continue
        for ln in ax.get_lines():
            if not ln.get_visible() or ln.get_transform() != ax.transData:
                continue
            xy = ln.get_xydata()
            if len(xy) == 0:
                continue
            disp = ax.transData.transform(xy)
            bb = ax.bbox
            for p in disp:
                if bb.x0 - 1 <= p[0] <= bb.x1 + 1 and bb.y0 - 1 <= p[1] <= bb.y1 + 1:
                    pts.append(p)
        for col in ax.collections:
            offs = col.get_offsets()
            if col.get_visible() and len(offs) and col.get_offset_transform() == ax.transData:
                pts.extend(ax.transData.transform(offs))
    on_data = []
    for t, b in boxes:
        if id(t) in owner:
            continue
        for p in pts:
            if b.x0 + 1 < p[0] < b.x1 - 1 and b.y0 + 1 < p[1] < b.y1 - 1:
                on_data.append(t.get_text()[:40])
                break
    return {
        "text_on_data": on_data,
        "size_in": [round(fig.get_figwidth(), 2), round(fig.get_figheight(), 2)],
        "size_mm": [round(fig.get_figwidth() * 25.4), round(fig.get_figheight() * 25.4)],
        "n_text": len(texts),
        "min_font_pt": min(sizes),
        "max_font_pt": max(sizes),
        "clipped_text": clipped,
        "text_overlaps": overlaps,
    }


def main() -> None:
    out = {}
    for name, mod in MODULES.items():
        captured: dict = {}
        m = importlib.import_module(mod)
        m.save = lambda fig, n, _c=captured, **kw: _c.__setitem__("fig", fig)
        m.main()
        out[name] = inspect(captured["fig"])
        plt.close("all")
    (_common.OUT_DIR / "figure_qc_auto.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False) + "\n"
    )
    for k, v in out.items():
        print(
            k,
            v["size_mm"],
            "min",
            v["min_font_pt"],
            "clipped",
            len(v["clipped_text"]),
            "overlaps",
            len(v["text_overlaps"]),
            "on-data",
            len(v["text_on_data"]),
        )
        for o in v["text_on_data"]:
            print("    on data:", o)
        for o in v["text_overlaps"]:
            print("    overlap:", o)
        for c in v["clipped_text"]:
            print("    clipped:", c)


if __name__ == "__main__":
    main()
