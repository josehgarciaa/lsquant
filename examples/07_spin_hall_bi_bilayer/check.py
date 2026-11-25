#!/usr/bin/env python3
"""Check Bi topology and Hall outputs; see shared light-input checks.

Pitfall: the quick Jackson curve does not resolve the insulating plateau.
Production is gated against the exact oracle even if convergence is pending.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from heavy_checks import guarded, spin_example

if __name__ == '__main__':
    sys.exit(guarded(lambda: spin_example('bi')))
