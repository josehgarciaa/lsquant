"""Diagnostics for legacy partial basis averages labelled as exact traces."""

import json
import os
import warnings


def trace_coverage(states, dimension, *, receipt_path=None, warn=True):
    """Warn about a truncated local average and return its coverage receipt.

    Pitfall: the legacy engine normalises by states, not dimension. Keeping
    its labels and arithmetic preserves existing moment and curve files.
    """
    if states < dimension:
        if warn:
            warnings.warn(
                f"Exact trace requested: computing an average over the first {states} of "
                f"{dimension} basis functions, not the trace",
                RuntimeWarning, stacklevel=2,
            )
        coverage = dict(states=int(states), dimension=int(dimension))
        if receipt_path and os.path.isfile(receipt_path):
            with open(receipt_path, encoding="utf-8") as fh:
                resolved = json.load(fh)
            resolved["expansion"]["trace_coverage"] = coverage
            coverages = resolved.setdefault("receipt", {}).setdefault("trace_coverages", [])
            if coverage not in coverages:
                coverages.append(coverage)
            with open(receipt_path, "w", encoding="utf-8") as fh:
                json.dump(resolved, fh, indent=2)
                fh.write("\n")
        return coverage
    return None
