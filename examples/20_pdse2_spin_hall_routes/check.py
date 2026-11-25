#!/usr/bin/env python3
"""Check the four published PdSe2 routes using light tabulated inputs.

Pitfall: dense-mesh exact curves are not same-mesh KPM references.
Use --data-only to test the page independently of any new calculation.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from heavy_checks import guarded, routes_example

if __name__ == '__main__':
    sys.exit(guarded(lambda: routes_example(20)))
