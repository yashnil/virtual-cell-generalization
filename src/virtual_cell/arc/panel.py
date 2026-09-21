"""Mapping a predicted response onto the official Arc gene panel.

The challenge scores a fixed 18,533-gene axis. Any response we predict is
defined on whatever gene space the public source data supplies. Genes in the
panel therefore fall into categories that differ in *why* we have no value for
them, and those categories must stay visible: a gene we measured and predicted
to be unchanged is a claim, whereas a gene we never measured is an absence of
evidence. Collapsing both to zero throws that distinction away.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np
import pandas as pd

__all__ = ["Support", "PanelMap", "build_panel_map", "project_response"]


class Support(StrEnum):
    """Why a panel gene does or does not carry a predicted response."""

    #: Measured in a source context and carrying a transferable response.
    PREDICTED = "predicted"
    #: Measured in a source context, but the perturbation was never applied
    #: there, so no response is estimable for this perturbation.
    MEASURED_NO_RESPONSE = "measured_no_response"
    #: Absent from every source gene space. Nothing is known about it.
    UNMEASURED = "unmeasured"


@dataclass(frozen=True)
class PanelMap:
    """An alignment between a source gene space and the Arc panel."""

    panel_genes: pd.Index
    source_genes: pd.Index
    #: For each panel gene, its position in ``source_genes``, or -1.
    source_index: np.ndarray
    support: np.ndarray

    def __post_init__(self) -> None:
        n = len(self.panel_genes)
        if self.source_index.shape != (n,) or self.support.shape != (n,):
            raise ValueError("source_index and support must be one entry per panel gene")

    @property
    def n_panel(self) -> int:
        return len(self.panel_genes)

    def counts(self) -> pd.Series:
        """How many panel genes fall in each support category."""
        return (
            pd.Series(self.support).value_counts().reindex([s.value for s in Support], fill_value=0)
        )

    def mask(self, support: Support) -> np.ndarray:
        return self.support == support.value


def build_panel_map(
    panel_genes: pd.Index | list[str],
    source_genes: pd.Index | list[str],
    *,
    responsive: np.ndarray | None = None,
) -> PanelMap:
    """Align ``source_genes`` to ``panel_genes`` and categorise every panel gene.

    ``responsive`` optionally marks, per *source* gene, whether a response is
    actually estimable (e.g. the perturbation was applied in that context). A
    source gene that is present but not responsive lands in
    :attr:`Support.MEASURED_NO_RESPONSE` rather than :attr:`Support.PREDICTED`.
    """
    panel = pd.Index(panel_genes).astype(object)
    source = pd.Index(source_genes).astype(object)
    if panel.has_duplicates:
        raise ValueError("panel gene names must be unique")

    if source.has_duplicates:
        # Keep the first occurrence; record nothing silently.
        source = source[~source.duplicated()]

    pos = pd.Series(np.arange(len(source)), index=source)
    idx = pos.reindex(panel).to_numpy()
    source_index = np.where(np.isnan(idx), -1, np.nan_to_num(idx, nan=-1.0)).astype(np.int64)

    support = np.full(len(panel), Support.UNMEASURED.value, dtype=object)
    present = source_index >= 0
    if responsive is None:
        support[present] = Support.PREDICTED.value
    else:
        resp = np.asarray(responsive, dtype=bool)
        if resp.shape != (len(source),):
            raise ValueError("responsive must have one entry per (deduplicated) source gene")
        hit = resp[source_index[present]]
        where = np.flatnonzero(present)
        support[where[hit]] = Support.PREDICTED.value
        support[where[~hit]] = Support.MEASURED_NO_RESPONSE.value

    return PanelMap(
        panel_genes=panel, source_genes=source, source_index=source_index, support=support
    )


def project_response(
    response: np.ndarray,
    panel_map: PanelMap,
    *,
    fill: float = 0.0,
) -> np.ndarray:
    """Lift a source-space response vector onto the panel axis.

    ``fill`` is used for every panel gene without a predicted value. It is an
    explicit argument with no default hiding: choosing 0.0 is a *decision* that
    unsupported genes are predicted unchanged, and the caller states it.
    """
    response = np.asarray(response, dtype=np.float64)
    if response.shape[-1] != len(panel_map.source_genes):
        raise ValueError(
            f"response has {response.shape[-1]} genes, expected {len(panel_map.source_genes)}"
        )
    out = np.full(response.shape[:-1] + (panel_map.n_panel,), fill, dtype=np.float64)
    take = panel_map.mask(Support.PREDICTED)
    out[..., take] = response[..., panel_map.source_index[take]]
    return out
