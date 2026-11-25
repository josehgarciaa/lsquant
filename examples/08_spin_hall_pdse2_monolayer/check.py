#!/usr/bin/env python3
"""Check PdSe2 topology and Hall outputs; see shared light-input checks.

The quick run file uses a 12x12 mesh; the calibration is for that mesh.
The upper-gap production value cannot calibrate conductivity at E_F = 0.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from heavy_checks import guarded, spin_example

if __name__ == '__main__':
    sys.exit(guarded(lambda: spin_example('pdse2')))
