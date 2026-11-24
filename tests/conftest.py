"""Thread count of the unit suite.

Problem. ``import lsquant`` sets ``OMP_NUM_THREADS`` to every physical core the process may
use when the variable is unset. The unit tests work on matrices of a few thousand rows, and
the engine opens one OpenMP parallel region per sparse product and per vector update. On a
node that other jobs share, a team as large as the machine then spends its time waiting for
cores: ``test_golden_chain_plans_are_identical_to_g1_g2`` (chain of 2048 sites, 64 moments)
takes 14 s with 32 threads and more than two minutes, which looks like a hang, with 64
threads on a 64-core node at load 35 to 60.

Approach. Before any test module imports the package, default ``OMP_NUM_THREADS`` to at most
four threads. A value set by the caller always wins, so a thread-scaling check still runs
with ``OMP_NUM_THREADS=64 python3 -m pytest tests/unit``.

Assumption. The moments do not depend on the thread count: the gate compares the same
goldens byte for byte at 1, 8, 16 and 64 threads (``scripts/check_goldens.sh``).

Pitfall. This file must not import ``lsquant``: the variable is read when libgomp loads.
"""
from __future__ import annotations

import os

#: Largest default team of the unit suite; the tests are too small to use more.
MAX_DEFAULT_THREADS = 4


def default_thread_count(cpu_count=None) -> int:
    """Return the default ``OMP_NUM_THREADS`` of the suite: 1 to ``MAX_DEFAULT_THREADS``."""
    if cpu_count is None:
        try:
            cpu_count = len(os.sched_getaffinity(0))
        except AttributeError:  # platforms without CPU affinity
            cpu_count = os.cpu_count() or 1
    return max(1, min(MAX_DEFAULT_THREADS, int(cpu_count)))


os.environ.setdefault("OMP_NUM_THREADS", str(default_thread_count()))
