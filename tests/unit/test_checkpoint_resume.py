"""Resume a genuine partial trace without changing the calculation identity.

The diagnostics hook throws only after the engine writes a completed state's
partial tensor. The saved state count proves that the second invocation skips
finished states rather than starting over. All input and output live in tmp_path;
the three-site, eight-moment calculation needs only constant-sized tensors.

Pitfall: reducing expansion.states creates a different fingerprint and cannot
represent an interruption of this calculation. Diagnostics are process-global,
so restore them even when an assertion fails. Run with OMP_NUM_THREADS=1.
"""
from __future__ import annotations

from pathlib import Path
import struct
import sys

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python"))
pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant import _core  # noqa: E402
from lsquant.runfile import execute_runfile  # noqa: E402


@pytest.mark.parametrize("finished", [1, 2])
def test_partial_checkpoint_resumes_byte_identically(tmp_path, finished):
    """Compare resumed and fresh moment bytes with the same seed and four states.

    Read only the fixed checkpoint header, using the writer's native byte order
    and explicit-width fields. A second interruption at k+1 proves forward
    progress; the final run must remove the completed checkpoint. No checkpoint
    bytes are fabricated or edited, and the calculation fingerprint stays fixed.
    """
    (tmp_path / "toy.HAM.CSR").write_text(
        "3 7\n0.2 0  0.4 0  0.4 0  -0.3 0  0.5 0  0.5 0  0.1 0\n"
        "0 1 0 1 2 1 2\n0 2 5 7\n"
    )
    raw = {
        "model": {"kind": "csr", "dir": str(tmp_path), "label": "toy", "operators": []},
        "basis": {"name": "real"},
        "word": "Tr[delta(E-H)]",
        "expansion": {"moments": 8, "window": [-2.0, 2.0], "trace": "random",
                      "states": 4, "seed": 17, "checkpoint": True, "state_block": 1},
        "output": {"legacy_layout": True, "write_chebmom": False},
    }
    resumed_dir = tmp_path / "resumed"
    header = struct.Struct("=8sIIIqqq")
    try:
        checkpoint = None
        for count in (finished, finished + 1):
            _core.ConfigureEngineDiagnostics(False, False, count)
            with pytest.raises(RuntimeError, match=f"after {count} states"):
                execute_runfile(raw, out_dir_override=str(resumed_dir), repo_root=str(REPO))
            paths = list(resumed_dir.glob("checkpoint_v2_*.lsqckpt"))
            assert len(paths) == 1
            if checkpoint is not None:
                assert paths[0] == checkpoint
            checkpoint = paths[0]
            magic, version, rank, ndims, elements, total, done = header.unpack_from(
                checkpoint.read_bytes()
            )
            assert (magic, version, rank, ndims, elements) == (b"LSQCKPT\n", 1, 1, 1, 8)
            assert total == 4
            assert done == count < total
        _core.ConfigureEngineDiagnostics(False, False, -1)
        resumed = execute_runfile(raw, out_dir_override=str(resumed_dir), repo_root=str(REPO))
        fresh = execute_runfile(raw, out_dir_override=str(tmp_path / "fresh"), repo_root=str(REPO))
        resumed_bytes = [Path(p).read_bytes() for p in resumed["moment_files"]]
        fresh_bytes = [Path(p).read_bytes() for p in fresh["moment_files"]]
        assert len(resumed_bytes) == len(fresh_bytes) == 1
        assert resumed_bytes == fresh_bytes
        assert not checkpoint.exists()
    finally:
        _core.ConfigureEngineDiagnostics(False, False, -1)
