"""Observe binding in concurrent CLI processes, with two OpenMP threads each.

Problem. Decision-function tests cannot detect binding applied when the native
runtime loads. A temporary startup module sets affinity before the real
``python -m lsquant`` import and registers an exit probe. The probe calls the
GNU parallel-region ABI through the dependency of the loaded extension, so it
observes actual OpenMP participants rather than unrelated Python/BLAS threads.
All processes rendezvous inside their regions using pipes, without sleeps or
speed assertions. Work and storage are O(processes * threads).

Pitfall. Importing the package before setting affinity or binding variables
would test an already initialized runtime. No package import occurs here.
Warning. The GNU ABI is optional; other OpenMP implementations are skipped
because the extension has no public per-participant affinity callback.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import selectors
import subprocess
import sys
import time

import pytest

pytestmark = pytest.mark.skipif(
    sys.platform != "linux" or not hasattr(os, "sched_setaffinity"),
    reason="Linux thread affinity is required",
)

REPO = Path(__file__).resolve().parents[2]
PREFIX = "AFFINITY_PROBE "

# This module is imported by Python before lsquant.__init__, without replacing
# configure_threads, CLI detection, or the extension. The words command avoids
# numerical work and leaves its normal import/entry-point behavior intact.
STARTUP = r'''
import atexit
import ctypes
import io
import json
import os
import sys
import threading
import traceback

# Keep normal CLI text out of the pipe protocol. In particular, text buffering
# must not hide a ready record behind output already prefetched by readline.
protocol = sys.stdout
sys.stdout = io.StringIO()
mask = json.loads(os.environ["AFFINITY_TEST_MASK"])
if mask is not None:
    os.sched_setaffinity(0, mask)
initial = sorted(os.sched_getaffinity(0))
owned = len(initial) < (os.cpu_count() or len(initial))

def emit(value):
    print("AFFINITY_PROBE " + json.dumps(value), file=protocol, flush=True)

def probe():
    try:
        core = sys.modules.get("lsquant._core")
        if core is None:
            raise RuntimeError("CLI did not load lsquant._core")
        lib = ctypes.CDLL(core.__file__)
        try:
            parallel = lib.GOMP_parallel
            thread_num = lib.omp_get_thread_num
            team_size = lib.omp_get_num_threads
            proc_bind = lib.omp_get_proc_bind
        except AttributeError:
            emit({"skip": "Loaded extension lacks the GNU parallel-region ABI; "
                           "no public per-thread affinity callback is exposed"})
            return
        thread_num.restype = ctypes.c_int
        team_size.restype = ctypes.c_int
        proc_bind.restype = ctypes.c_int
        callback_type = ctypes.CFUNCTYPE(None, ctypes.c_void_p)
        parallel.argtypes = [callback_type, ctypes.c_void_p,
                             ctypes.c_uint, ctypes.c_uint]
        parallel.restype = None
        gate = threading.Barrier(2)
        samples = []
        errors = []

        @callback_type
        def participant(_):
            try:
                size = team_size()
                if size != 2:
                    raise RuntimeError("Expected two OpenMP threads, got " + str(size))
                tid = threading.get_native_id()
                samples.append({"index": thread_num(), "tid": tid,
                                "cpus": sorted(os.sched_getaffinity(tid))})
                gate.wait(timeout=4)
                if thread_num() == 0:
                    emit({"ready": True})
                    if sys.stdin.readline() != "go\n":
                        raise RuntimeError("Missing parallel-region release")
                gate.wait(timeout=4)
            except BaseException:
                errors.append(traceback.format_exc())
                gate.abort()

        parallel(participant, None, 0, 0)
        if errors:
            raise RuntimeError("\n".join(errors))
        emit({"initial": initial, "owned": owned, "threads": samples,
              "active_binding": proc_bind(),
              "bind": os.environ.get("OMP_PROC_BIND"),
              "places": os.environ.get("OMP_PLACES")})
    except BaseException:
        traceback.print_exc()
        os._exit(1)

atexit.register(probe)
'''


def _concurrent_affinities(tmp_path, masks, choice):
    """Keep all CLI regions alive together, then collect their participants.

    A shared deadline bounds failed synchronization, rather than asserting a
    runtime speed. Every child is reaped even after startup errors or skips.
    Pitfall: inherited placement settings override application defaults, so
    remove them, including vendor affinity controls, before each fresh import.
    """
    (tmp_path / "sitecustomize.py").write_text(STARTUP)
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(("OMP_", "GOMP_", "KMP_"))
           and key != "LSQUANT_BIND_THREADS"}
    env.update(OMP_NUM_THREADS="2", OMP_DYNAMIC="FALSE", OMP_THREAD_LIMIT="2",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1")
    # Preserve the configured extension search path, including relative entries,
    # even though children create any incidental files only under tmp_path.
    paths = [str(tmp_path), str(REPO / "python")]
    paths.extend(str(Path(p or ".").resolve()) for p in
                 os.environ.get("PYTHONPATH", "").split(os.pathsep))
    env["PYTHONPATH"] = os.pathsep.join(paths)
    if choice is not None:
        env["LSQUANT_BIND_THREADS"] = choice
    children = []
    results = []
    deadline = time.monotonic() + 4.5
    try:
        with selectors.DefaultSelector() as selector:
            for index, mask in enumerate(masks):
                err = (tmp_path / f"stderr-{index}.txt").open("w+")
                children.append((None, err))
                child = subprocess.Popen(
                    [sys.executable, "-m", "lsquant", "words", "list"],
                    cwd=tmp_path,
                    env=dict(env, AFFINITY_TEST_MASK=json.dumps(mask)),
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err,
                    text=True, bufsize=1,
                )
                children[-1] = (child, err)
                selector.register(child.stdout, selectors.EVENT_READ, child)
            ready = set()
            while len(ready) < len(children):
                events = selector.select(max(0, deadline - time.monotonic()))
                assert events, "CLI parallel-region rendezvous did not complete"
                for key, _ in events:
                    line = key.fileobj.readline()
                    if not line:
                        _, err = next(pair for pair in children if pair[0] is key.data)
                        err.seek(0)
                        pytest.fail("CLI exited before rendezvous: " + err.read())
                    if line.startswith(PREFIX):
                        message = json.loads(line[len(PREFIX):])
                        if "skip" in message:
                            pytest.skip(message["skip"])
                        assert message == {"ready": True}, message
                        ready.add(key.data.pid)
                        selector.unregister(key.fileobj)
            for child, _ in children:
                child.stdin.write("go\n")
                child.stdin.flush()
            for child, err in children:
                output, _ = child.communicate(timeout=max(0.01, deadline - time.monotonic()))
                err.seek(0)
                assert child.returncode == 0, err.read()
                records = [json.loads(line[len(PREFIX):])
                           for line in output.splitlines() if line.startswith(PREFIX)]
                assert len(records) == 1, output
                result = records[0]
                assert len(result["threads"]) == 2, result
                assert {t["index"] for t in result["threads"]} == {0, 1}, result
                assert len({t["tid"] for t in result["threads"]}) == 2, result
                results.append(result)
        return results
    finally:
        for child, err in children:
            if child is not None:
                if child.poll() is None:
                    child.kill()
                child.wait()
                child.stdin.close()
                child.stdout.close()
            err.close()


@pytest.mark.parametrize("count", [2, 4])
@pytest.mark.parametrize("choice", [None, "0"], ids=["automatic", "disabled"])
def test_shared_cli_threads_are_not_pinned(tmp_path, count, choice):
    """Shared commands leave every actual participant free to use multiple CPUs."""
    allowed = os.sched_getaffinity(0)
    if len(allowed) < 2:
        pytest.skip("At least two allowed CPUs are needed to detect single-CPU binding")
    if choice is None and len(allowed) < (os.cpu_count() or len(allowed)):
        pytest.skip("owns_its_cores() holds for the inherited affinity; automatic binding is valid")
    for result in _concurrent_affinities(tmp_path, [None] * count, choice):
        if choice is None:
            assert not result["owned"], result
        assert result["bind"] is None, result
        assert result["active_binding"] == 0, result
        for thread in result["threads"]:
            assert set(thread["cpus"]) == set(result["initial"]), result
            assert len(thread["cpus"]) > 1, result


def test_forced_cli_threads_respect_disjoint_masks(tmp_path):
    """Explicit binding stays within each process's two-CPU allocation."""
    allowed = sorted(os.sched_getaffinity(0))
    if len(allowed) < 4:
        pytest.skip("Disjoint two-CPU allocations require at least four allowed CPUs")
    masks = [allowed[:2], allowed[2:4]]
    results = _concurrent_affinities(tmp_path, masks, "1")
    used = []
    for mask, result in zip(masks, results):
        assert result["initial"] == mask, result
        assert result["bind"] == "spread", result
        assert result["places"] == "cores", result
        assert result["active_binding"] != 0, result
        cpus = set()
        for thread in result["threads"]:
            assert thread["cpus"], result
            assert set(thread["cpus"]) <= set(mask), result
            cpus.update(thread["cpus"])
        used.append(cpus)
    assert used[0].isdisjoint(used[1]), results
