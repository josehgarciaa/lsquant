"""CPU affinity and kernel tuning defaults applied before loading OpenMP or NumPy."""
from __future__ import annotations

import os
import sys
from pathlib import Path

DEFAULT_SPMV_MIN_NNZ = "16384"
DEFAULT_KERNEL_BATCH = "0"


def physical_core_count() -> int:
    """Count distinct physical cores in the process CPU affinity, including cpusets."""
    try:
        allowed = os.sched_getaffinity(0)
    except AttributeError:
        return os.cpu_count() or 1
    cores = set()
    for cpu in allowed:
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        try:
            cores.add(((topology / "physical_package_id").read_text().strip(),
                       (topology / "core_id").read_text().strip()))
        except OSError:
            cores.add(("cpu", str(cpu)))
    return max(1, len(cores))


def launched_as_cli(argv=None) -> bool:
    """True when this process is the ``lsquant`` command (console script or ``python -m lsquant``).

    ``sys.orig_argv`` keeps ``-m lsquant`` while the package is imported; ``sys.argv[0]`` is
    still ``-m`` at that point.
    """
    args = list(argv if argv is not None else getattr(sys, "orig_argv", sys.argv))
    for flag, value in zip(args, args[1:]):
        if flag == "-m" and value in ("lsquant", "lsquant.cli"):
            return True
    return bool(args) and any(Path(a).name == "lsquant" for a in args[:2])


def configure_threads(argv=None) -> None:
    """Thread count for every process; core binding only for the ``lsquant`` command.

    With ``OMP_PROC_BIND`` set, libgomp binds the initial thread to the first place when it
    loads, so every Python process that imports the package (scripts, notebooks,
    multiprocessing workers) would run on one core. Spread binding, which keeps engine threads
    off a single CCD, is therefore applied only when the process is the command itself. Values
    set by the caller always win.
    """
    if launched_as_cli(argv):
        os.environ.setdefault("OMP_PROC_BIND", "spread")
        os.environ.setdefault("OMP_PLACES", "cores")
    if "OMP_NUM_THREADS" not in os.environ:
        os.environ["OMP_NUM_THREADS"] = str(physical_core_count())


def configure_kernels(core) -> None:
    """Send numeric tuning and diagnostic settings once; the engine never reads the environment.

    Diagnostics (both profiles): ``LSQ_ENGINE_MEM`` and ``LSQ_ENGINE_TIMING`` (set to any value)
    print the allocation and the phase timings of the native rank-2 path to stderr;
    ``LSQ_CHECKPOINT_ABORT_AFTER=k`` stops a checkpointed run right after the k-th checkpoint
    write, for the resume test.
    """
    if hasattr(core, "ConfigureEngineDiagnostics"):
        core.ConfigureEngineDiagnostics(
            "LSQ_ENGINE_MEM" in os.environ,
            "LSQ_ENGINE_TIMING" in os.environ,
            int(os.environ.get("LSQ_CHECKPOINT_ABORT_AFTER", "-1")),
        )
    if core.build_profile == "native":
        core.ConfigurePerformance(
            int(os.environ.get("LSQUANT_SPMV_MIN_NNZ", DEFAULT_SPMV_MIN_NNZ)),
            int(os.environ.get("LSQUANT_KERNEL_BATCH", DEFAULT_KERNEL_BATCH)),
        )
