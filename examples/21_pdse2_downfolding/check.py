#!/usr/bin/env python3
"""Check published downfolding plateaus and supplied SIESTA output curves.

Pitfall: reduced and original exact curves have different occupied states.
The shared implementation compares each curve with its own reference.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from heavy_checks import guarded, routes_example

if __name__ == '__main__':
    sys.exit(guarded(lambda: routes_example(21)))
