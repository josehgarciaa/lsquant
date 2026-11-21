#!/usr/bin/env python3
"""Drive the LSQuant MCP server over stdio with the ``mcp`` client library.

The script starts ``python -m lsquant.mcp_server --root <repo>`` as a child
process, opens an MCP client session on its stdin/stdout and calls the five
tools the server exposes:

    words_describe   parse and compile an observable word
    run              validate a run file and start it in the background
    status           poll the background run
    inspect          read the run directory (resolved run file + moment metadata)
    reconstruct      rebuild curves from the stored moments at a new broadening

The calculation is the Bloch-basis graphene run of example 10 (about a second
of compute). All paths the server touches must lie inside --root; the output
goes to <root>/out/mcp_session.

Requires: pip install mcp (tested with mcp 2.3.0).

Usage:
    python examples/18_mcp_session/mcp_session.py
"""
import asyncio
import json
import os
import sys
import time

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(_HERE))
OUT = "out/mcp_session"


def is_error(result):
    """Tool-error flag (``is_error`` in mcp 2.x, ``isError`` in mcp 1.x)."""
    return bool(getattr(result, "is_error", getattr(result, "isError", False)))


def text(result):
    """Return the text payload of a tool result, raising on a tool error."""
    payload = "".join(getattr(c, "text", "") for c in result.content)
    if is_error(result):
        raise RuntimeError(payload)
    return payload


async def main():
    with open(os.path.join(ROOT, "examples", "10_bases_compare", "run_bloch.json")) as fh:
        run = json.load(fh)
    run["output"]["dir"] = OUT
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [os.path.join(ROOT, "python")] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    server = StdioServerParameters(command=sys.executable,
                                   args=["-m", "lsquant.mcp_server", "--root", ROOT], env=env)
    async with stdio_client(server) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            tools = await session.list_tools()
            print("tools:", sorted(t.name for t in tools.tools))

            info = json.loads(text(await session.call_tool(
                "words_describe", {"word": "Tr[VX delta(E-H) VY dG-(E)/dE]"})))
            print("words_describe: shape =", info["reconstruction"]["shape_id"],
                  "operators =", info["operators"], "formulas =", info["reconstruction"]["formulas"],
                  "plan =", len(info["plan"]["instructions"]), "instructions")

            started = json.loads(text(await session.call_tool("run", {"runfile_json": json.dumps(run)})))
            print("run:", started)
            t0 = time.time()
            while True:
                st = json.loads(text(await session.call_tool("status", {"run_id": started["run_id"]})))
                if st["status"] != "running":
                    break
                await asyncio.sleep(1.0)
            print(f"status: {st['status']} (returncode {st['returncode']}) after {time.time() - t0:.1f} s")
            if st["status"] != "succeeded":
                print(open(st["log"]).read())
                return 1

            ins = json.loads(text(await session.call_tool("inspect", {"path": OUT})))
            print("inspect:", ins["kind"], "basis =", ins["resolved_runfile"]["basis"]["name"],
                  "moment files =", [(os.path.basename(m["path"]), m["dims"]) for m in ins["moment_files"]])

            dos_file = OUT + "/SpectralOp1graphene_bloch_KPM_M128_stateexact.chebmom1D"
            rec = json.loads(text(await session.call_tool(
                "reconstruct", {"output_dir": dos_file, "formula": "dos", "broadening_meV": 300.0})))
            print("reconstruct:", rec["formula"], rec["broadening_meV"], "meV ->",
                  [os.path.basename(p) for p in rec["curve_files"]])

            res = await session.call_tool("inspect", {"path": "/etc"})
            print("root guard: inspect('/etc') refused =", is_error(res), "|", text_or_error(res)[:100])
            return 0 if is_error(res) else 1


def text_or_error(result):
    return "".join(getattr(c, "text", "") for c in result.content)


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main()))
