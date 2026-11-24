"""FastMCP stdio server exposing LSQUANT word compilation, run dispatch, reconstruction, and inspection."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import os
import subprocess
import sys
from typing import Any
import uuid

from .runfile import (
    RunFile,
    inspect_target,
    reconstruct_target,
    resolve_generator_path,
)
from .words import describe


def _default_root() -> str:
    env_root = os.environ.get("LSQUANT_MCP_ROOT")
    if env_root:
        return os.path.realpath(os.path.abspath(env_root))
    return os.path.realpath(
        os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    )


_CONFIGURED_ROOT: str = _default_root()


def set_configured_root(root_dir: str) -> str:
    """Set and return the canonical root directory for MCP file access."""
    global _CONFIGURED_ROOT  # noqa: PLW0603
    _CONFIGURED_ROOT = os.path.realpath(os.path.abspath(root_dir))
    return _CONFIGURED_ROOT


def get_configured_root() -> str:
    """Return the canonical root directory for MCP file access."""
    return _CONFIGURED_ROOT


def _resolve_within_root(
    path_str: str, root: str | None = None, *, key: str = "path",
    supplied: str | None = None,
) -> str:
    """Resolve ``path_str`` and refuse access if it lies outside ``root``."""
    base_root = os.path.realpath(root if root is not None else _CONFIGURED_ROOT)
    candidate = (
        path_str
        if os.path.isabs(path_str)
        else os.path.join(base_root, path_str)
    )
    resolved = os.path.realpath(os.path.abspath(candidate))
    try:
        common = os.path.commonpath([base_root, resolved])
    except ValueError:
        common = None
    if common != base_root:
        raise ValueError(
            f"Access denied: {key} '{path_str if supplied is None else supplied}' "
            "is outside the server root"
        )
    return resolved


# All explicit run-file paths belong here, including connector parameters even
# when supplied on another model kind. The mode records the executor's base.
_RUN_PATH_FIELDS = (
    ("model.dir", "resolved"),
    ("model.csr_dir", "resolved"),
    ("model.file", "resolved"),
    ("model.hr_file", "resolved"),
    ("model.positions", "resolved"),
    ("model.atom_map", "resolved"),
    ("model.name", "generator"),
    ("model.params.overlap_file", "resolved"),
    ("model.params.stru_file", "resolved"),
    ("model.params.orbital_dir", "resolved"),
    ("model.params.position_file", "resolved"),
    ("model.params.dipole_cache", "cwd"),
    ("model.params.hamiltonian_file", "model_dir"),
    ("output.dir", "cwd"),
)


def _guard_runfile_paths(rf: RunFile, raw: dict[str, Any], root: str) -> str:
    """Check explicit inputs and outputs before creating files or dispatching.

    RunFile resolves ordinary inputs against the root, then the source file's
    directory. Connector exceptions follow the subprocess working directory or
    model.dir instead. One field list keeps these checks independent of connector
    dispatch; string positions/atom maps are paths, while inline arrays are not.

    Pitfall: report only the caller's spelling, never a symlink's private target.
    The root itself is allowed, preserving the default output.dir of ".".
    Warning: this is a path check, not protection against concurrent filesystem
    changes between validation and execution.
    """
    output_dir = root
    for key, mode in _RUN_PATH_FIELDS:
        parts = key.split(".")
        value: Any = getattr(rf, parts[0])
        original: Any = raw.get(parts[0], {})
        for part in parts[1:]:
            value = value.get(part) if isinstance(value, dict) else None
            original = original.get(part) if isinstance(original, dict) else None
        if not isinstance(value, str):
            continue
        candidate = value
        if mode == "generator":
            candidate = resolve_generator_path(value, root)
        elif mode == "model_dir":
            candidate = os.path.join(rf.model.get("dir") or root, value)
        resolved = _resolve_within_root(
            candidate, root, key=key, supplied=str(original) if original is not None else value,
        )
        if key == "model.params.dipole_cache" and not value.endswith(".npz"):
            # numpy.savez appends this suffix when writing a cache miss.
            _resolve_within_root(candidate + ".npz", root, key=key, supplied=value)
        if key == "output.dir":
            output_dir = resolved
    return output_dir


@dataclass
class _DispatchedRun:
    run_id: str
    process: subprocess.Popen[bytes]
    log_path: str
    output_dir: str
    log_handle: Any


_ACTIVE_RUNS: dict[str, _DispatchedRun] = {}


def tool_words_describe(word: str) -> str:
    """Validate and describe an observable word, returning JSON identical to ``lsquant words describe``."""
    info = describe(word)
    return json.dumps(info, indent=2)


def tool_run(runfile_json: str | dict[str, Any]) -> str:
    """Validate a run file and dispatch ``python -m lsquant run`` in a non-blocking subprocess."""
    root = get_configured_root()
    if isinstance(runfile_json, str):
        stripped = runfile_json.strip()
        if stripped.startswith("{"):
            raw = json.loads(stripped)
            source_path = None
        else:
            resolved_rf_path = _resolve_within_root(stripped, root, key="runfile_json")
            with open(resolved_rf_path, encoding="utf-8") as fh:
                raw = json.load(fh)
            source_path = resolved_rf_path
    elif isinstance(runfile_json, dict):
        raw = dict(runfile_json)
        source_path = None
    else:
        raise ValueError("runfile_json must be a JSON string, file path, or dict")

    rf = RunFile.from_dict(raw, source_path=source_path, repo_root=root)

    out_dir = _guard_runfile_paths(rf, raw, root)
    os.makedirs(out_dir, exist_ok=True)

    run_id = uuid.uuid4().hex[:12]
    runfile_copy = os.path.join(out_dir, f".mcp_run_{run_id}.json")
    log_path = os.path.join(out_dir, f"run_{run_id}.log")

    resolved_payload = rf.to_resolved_dict()
    resolved_payload["output"] = dict(resolved_payload["output"], dir=out_dir)
    with open(runfile_copy, "w", encoding="utf-8") as fh:
        json.dump(resolved_payload, fh, indent=2)
        fh.write("\n")

    log_fh = open(log_path, "wb")  # noqa: SIM115
    proc = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "lsquant",
            "run",
            runfile_copy,
            "--out",
            out_dir,
        ],
        cwd=root,
        stdout=log_fh,
        stderr=subprocess.STDOUT,
        env=dict(os.environ),
    )
    _ACTIVE_RUNS[run_id] = _DispatchedRun(
        run_id=run_id,
        process=proc,
        log_path=log_path,
        output_dir=out_dir,
        log_handle=log_fh,
    )
    return json.dumps(
        {
            "run_id": run_id,
            "log": log_path,
            "output_dir": out_dir,
        },
        indent=2,
    )


def tool_status(run_id: str) -> str:
    """Return the status of a previously dispatched background run."""
    if run_id not in _ACTIVE_RUNS:
        raise ValueError(f"Unknown run_id '{run_id}'")
    entry = _ACTIVE_RUNS[run_id]
    rc = entry.process.poll()
    if rc is None:
        state = "running"
    else:
        try:
            entry.log_handle.flush()
            entry.log_handle.close()
        except Exception:  # noqa: BLE001
            pass
        state = "succeeded" if rc == 0 else "failed"
    return json.dumps(
        {
            "run_id": entry.run_id,
            "status": state,
            "returncode": rc,
            "log": entry.log_path,
            "output_dir": entry.output_dir,
        },
        indent=2,
    )


def tool_reconstruct(
    output_dir: str,
    formula: str | None = None,
    broadening_meV: float | None = None,
    fermi_energies: list[float] | None = None,
) -> str:
    """Reconstruct curve(s) from an output directory or moment file within the configured root."""
    resolved_target = _resolve_within_root(output_dir, key="output_dir")
    res = reconstruct_target(
        resolved_target,
        formula=formula,
        broadening_meV=broadening_meV,
        fermi_energies=fermi_energies,
    )
    return json.dumps(res, indent=2)


def tool_inspect(path: str) -> str:
    """Inspect a moment file, run file, or directory within the configured root, returning JSON identical to ``lsquant inspect``."""
    resolved_path = _resolve_within_root(path, key="path")
    info = inspect_target(resolved_path)
    return json.dumps(info, indent=2)


def create_mcp_server(root_dir: str | None = None) -> Any:
    """Create and configure the FastMCP / MCPServer instance."""
    try:
        from mcp.server.mcpserver import MCPServer as _ServerClass
        from mcp.server.mcpserver.exceptions import ToolError as _ToolError
    except ImportError:
        from mcp.server.fastmcp import FastMCP as _ServerClass
        from mcp.server.fastmcp.exceptions import ToolError as _ToolError

    if root_dir is not None:
        set_configured_root(root_dir)

    mcp = _ServerClass("lsquant")

    @mcp.tool(name="words_describe")
    def words_describe(word: str) -> str:
        """Parse and compile an observable word, returning its AST, bytecode plan, operators, and reconstruction spec as JSON."""
        try:
            return tool_words_describe(word)
        except Exception as exc:  # noqa: BLE001
            raise _ToolError(str(exc)) from exc

    @mcp.tool(name="run")
    def run(runfile_json: str) -> str:
        """Validate a run file JSON string or path and launch the calculation in a background subprocess."""
        try:
            return tool_run(runfile_json)
        except Exception as exc:  # noqa: BLE001
            raise _ToolError(str(exc)) from exc

    @mcp.tool(name="status")
    def status(run_id: str) -> str:
        """Check the execution status of a background calculation started by ``run``."""
        try:
            return tool_status(run_id)
        except Exception as exc:  # noqa: BLE001
            raise _ToolError(str(exc)) from exc

    @mcp.tool(name="reconstruct")
    def reconstruct(
        output_dir: str,
        formula: str | None = None,
        broadening_meV: float | None = None,
        fermi_energies: list[float] | None = None,
    ) -> str:
        """Reconstruct spectral or transport curves from moments in ``output_dir``."""
        try:
            return tool_reconstruct(
                output_dir=output_dir,
                formula=formula,
                broadening_meV=broadening_meV,
                fermi_energies=fermi_energies,
            )
        except Exception as exc:  # noqa: BLE001
            raise _ToolError(str(exc)) from exc

    @mcp.tool(name="inspect")
    def inspect(path: str) -> str:
        """Inspect a Chebyshev moment file, run file, or output directory."""
        try:
            return tool_inspect(path)
        except Exception as exc:  # noqa: BLE001
            raise _ToolError(str(exc)) from exc

    return mcp


def main(argv: list[str] | None = None) -> int:
    """Entry point for ``python -m lsquant.mcp_server`` over stdio."""
    parser = argparse.ArgumentParser(prog="lsquant.mcp_server")
    parser.add_argument(
        "--root",
        default=None,
        help="Root directory restricting file access (defaults to LSQUANT_MCP_ROOT or repository root).",
    )
    args = parser.parse_args(argv)
    server = create_mcp_server(root_dir=args.root)
    server.run(transport="stdio")
    return 0


if __name__ == "__main__":
    sys.exit(main())
