# Example 18: driving LSQuant from an AI assistant through MCP

## The physical question

The Model Context Protocol (MCP) lets an AI assistant, such as Claude or any
other MCP client, call tools on a local machine. LSQuant ships an MCP server.
Through it, an assistant can check an observable word, launch a calculation,
watch it and read back the results. The server uses the same validation and
the same engine as the command line.

In this example a short Python script plays the role of the assistant. It
opens an MCP session with the server and calls every tool. The calculation is
the graphene DOS and $\sigma_{xx}$ of
[example 10](../10_bases_compare/README.md) in the Bloch basis. At the end we
change the broadening of the DOS from the stored moments and read the new
curve against the graphene bands.

## The server

The server is started with:

```bash
python3 -m lsquant.mcp_server --root /path/to/lsquant
```

It communicates over standard input and output and exposes five tools:

| Tool | Action | CLI equivalent |
|---|---|---|
| `words_describe(word)` | parses and compiles a word: shape, operators, formulas, plan | `lsquant words describe` |
| `run(runfile_json)` | validates a run file, given as JSON text or a path, and starts it in a background process; returns a `run_id` | `lsquant run` |
| `status(run_id)` | `running`, `succeeded` or `failed`, with the return code and log path | — |
| `inspect(path)` | the JSON of `lsquant inspect` | `lsquant inspect` |
| `reconstruct(output_dir, formula, broadening_meV, fermi_energies)` | rebuilds curves from stored moments | `lsquant reconstruct` |

All files that the server reads or writes must lie under its root. The root is
set by `--root`, or else by the `LSQUANT_MCP_ROOT` environment variable, and it
defaults to the repository. Relative paths are taken from the root.

An assistant registers the server with the command above. For Claude Code, the
registration reads:

```bash
claude mcp add lsquant -- python3 -m lsquant.mcp_server --root /path/to/lsquant
```

The server needs the `mcp` Python package, installed with `pip install mcp`.
It was tested with mcp 2.3.0.

## The session

[`mcp_session.py`](mcp_session.py) starts the server as a child process with
the `mcp` client library and opens a session:

```python
server = StdioServerParameters(command=sys.executable,
                               args=["-m", "lsquant.mcp_server", "--root", ROOT], env=env)
async with stdio_client(server) as (read, write):
    async with ClientSession(read, write) as session:
        await session.initialize()
```

We ask the server to describe the Hall word:

```python
await session.call_tool("words_describe", {"word": "Tr[VX delta(E-H) VY dG-(E)/dE]"})
```

We send the run file of example 10 as JSON text, with its output directory
moved to `out/mcp_session` under the root. We then poll `status` until the
background run ends:

```python
started = await session.call_tool("run", {"runfile_json": json.dumps(run)})
await session.call_tool("status", {"run_id": run_id})
```

We inspect the output directory, and we reconstruct the DOS from the stored
moments at 300 meV instead of 150 meV:

```python
await session.call_tool("inspect", {"path": "out/mcp_session"})
await session.call_tool("reconstruct", {"output_dir": dos_file, "formula": "dos",
                                        "broadening_meV": 300.0})
```

Finally we check that the server refuses a path outside its root. The whole
session runs with:

```bash
python3 examples/18_mcp_session/mcp_session.py
```

## The result

The script prints one line per tool call:

```text
tools: ['inspect', 'reconstruct', 'run', 'status', 'words_describe']
words_describe: shape = kubo_bastin operators = ['VX', 'VY'] formulas = ['bastin'] plan = 6 instructions
run: {'run_id': 'aa20198b933d', 'log': '<root>/out/mcp_session/run_aa20198b933d.log', 'output_dir': '<root>/out/mcp_session'}
status: succeeded (returncode 0) after 2.0 s
inspect: directory basis = bloch moment files = [('NonEqOpVX-VXgraphene_bloch_KPM_M128x128_stateexact.chebmom', [128, 128]), ...]
reconstruct: dos 300.0 meV -> ['meanSpectralOp1graphene_bloch_KPM_M128_stateexactJACKSON.dat']
root guard: inspect('/etc') refused = True | Error executing tool inspect: Access denied: path '/etc' resolves to '/etc' outside configured root
```

Every tool answers with the same JSON as its command-line equivalent. A tool
error comes back as an MCP error result, with `is_error` set to true and the
LSQuant message, and the session continues.

The `run` tool returns at once. The calculation runs in its own process and
writes its log to `run_<id>.log` in the output directory, and `status` reports
the outcome. The `reconstruct` tool, given a directory, applies one formula to
every moment file in it. Its default formula is `dos` for rank-1 files and
`bastin` for rank-2 files. A run with both a DOS word and a Kubo–Greenwood word
therefore needs the moment file and the formula given explicitly, as in this
session.

The reconstructed curve is the graphene DOS of example 10, now at 300 meV.
With the window $[-9, 9]$ eV, a 300 meV broadening uses
$\lceil \pi \cdot 9/0.3 \rceil = 95$ of the 128 stored moments. The 150 meV
request of the original run asked for more moments than were stored, so that
curve used all 128, a resolution of $\pi \cdot 9/128 \approx 0.22$ eV. The new
curve keeps the shape of the graphene bands: a DOS that vanishes linearly at
the Dirac point $E = 0$ and van Hove peaks at $E = \pm|t| = \pm 2.7$ eV. The
wider kernel lowers the van Hove peak from 90.8 to 71.9 states/eV, still at
$E = 2.70$ eV. It also fills the dip at the Dirac point from 0.03 to
1.3 states/eV. The number of states stays the same: both curves integrate to
512.0, the $2 \times 256$ states of the $16 \times 16$ supercell.

These values are for the whole supercell. The session prints no curve values:
its `reconstruct` line above gives only the name of the curve file. We
therefore measured the values by reconstructing the moments of example 10 at
300 meV on the command line.

## Runtime

On a shared 96-core node with one OpenMP thread, the whole session took about
5 s, including the 2 s background run.
