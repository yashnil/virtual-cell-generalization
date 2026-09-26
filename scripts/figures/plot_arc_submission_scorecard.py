"""Figure 10 — Arc submission v1 scorecard (ONLY after a real score exists).

Reads ``data/figure_sources/arc_submission_v1_scorecard.json``, which must be
filled from ``vcc status <entry> --json`` after scoring (schema:
``data/figure_sources/arc_submission_v1_scorecard.schema.json``). If that file
does not exist, the script writes the schema, draws **nothing**, and exits 0.
No placeholder or expected value is ever plotted.

Reference marks come from frozen artifacts only: 0 = mean-response baseline,
1 = replicate anchor, and the per-member control-emitting score from
``outputs/arc_bridge_v1/score_accounting.csv``.

Output: ``reports/figures/10_arc_submission_scorecard.{png,svg}``
"""

from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import pandas as pd

from virtual_cell.visualization import common, scorecard, style


def main() -> int:
    scorecard.write_schema()
    try:
        record = scorecard.load_scorecard()
    except scorecard.ScorecardError as exc:
        print(f"  no scorecard drawn: {exc}")
        return 0
    style.apply()
    acct = pd.read_csv(scorecard.ACCOUNTING_PATH).set_index("member")
    fig, (ax0, ax) = plt.subplots(
        2, 1, figsize=(style.FULL_WIDTH, 4.2), gridspec_kw={"height_ratios": [1, 3.2]}
    )
    control_overall = float(acct.control_submission_score.mean())
    ax0.barh(0, record["overall"], color=style.BLUE, edgecolor="black", height=0.5)
    ax0.plot(control_overall, 0, "|", color=style.VERMILION, ms=18, mew=2)
    ax0.set_yticks([0], ["Overall"])
    for a in (ax0, ax):
        style.reference_line(a, 0, axis="x")
        style.reference_line(a, 1, axis="x", color=style.GREY)
    for i, (key, _label, member) in enumerate(scorecard.MEMBERS):
        v = record[key]
        ax.barh(i, v, color=style.BLUE, edgecolor="black", height=0.6)
        ax.text(v, i, f" {v:+.3f}", va="center", fontsize=7)
        ax.plot(
            acct.loc[member, "control_submission_score"],
            i,
            "|",
            color=style.VERMILION,
            ms=14,
            mew=2,
        )
    ax.set_yticks(range(len(scorecard.MEMBERS)), [m[1] for m in scorecard.MEMBERS])
    ax.invert_yaxis()
    ax.set_xlabel("official normalized score (0 = mean-response baseline, 1 = replicate)")
    ax0.set_title(
        f"Arc submission v1 — partition {record['partition']}, panel {record['panel']}, "
        f"anchors {record['anchor_set']}",
        fontsize=9,
    )
    ax.plot(
        [],
        [],
        "|",
        color=style.VERMILION,
        ms=10,
        mew=2,
        label="control-emitting submission (frozen accounting)",
    )
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    common.save_figure(fig, "10_arc_submission_scorecard")
    print("  wrote reports/figures/10_arc_submission_scorecard.{png,svg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
