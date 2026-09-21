"""Biological priors that can represent a perturbation never measured.

The four-context study predicted a perturbation's conserved effect ``beta_p``
by *averaging measurements of that same perturbation* in other contexts. For a
perturbation that was never applied anywhere, that route does not exist, and
the only thing left is what is known about the gene itself independently of any
knockdown experiment.

This package assembles those representations and, just as importantly, records
for each one whether it could leak a perturbation outcome. A feature derived
from perturbation responses would make an unseen-perturbation benchmark
meaningless, so every source carries an explicit leakage classification in
:data:`~virtual_cell.priors.catalogue.CATALOGUE` and the audit that produced it
lives in ``data/provenance/``.
"""

from virtual_cell.priors import catalogue, features

__all__ = ["catalogue", "features"]
