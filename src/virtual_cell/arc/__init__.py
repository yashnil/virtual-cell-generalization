"""Arc Virtual Cell Challenge 2026 bridge.

Everything needed to turn a response prediction into a submittable object, and
to score one locally against a reference.

This package deliberately contains **no model**. It is the transport layer
between the frozen research findings and the challenge's output space:

* :mod:`~virtual_cell.arc.panel` maps a response defined on a source gene
  space onto the official 18,533-gene panel, with every gene assigned to an
  explicit support category rather than silently zero-filled.
* :mod:`~virtual_cell.arc.generate` turns a profile into raw integer counts.
* :mod:`~virtual_cell.arc.metrics` reimplements the six scored ``vcc2026``
  metrics so a candidate can be measured without submitting it.
"""

from virtual_cell.arc import generate, metrics, panel

__all__ = ["generate", "metrics", "panel"]
