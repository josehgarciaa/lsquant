# Connect an AI assistant

This guide shows you how to connect an AI assistant, such as Claude Code or
Claude Desktop, to LSQuant. The assistant communicates with LSQuant through
the **MCP server**, a [Model Context Protocol](https://modelcontextprotocol.io)
server that LSQuant ships as `python -m lsquant.mcp_server`. Through it, the
assistant can describe an observable word, start a run, follow it, inspect the
output and reconstruct curves, with the same validation as the command line.

The [MCP tool reference](../reference/mcp-tools.md) lists the five tools, their
arguments and their return values.
[Example 18](../../examples/18_mcp_session/README.md) runs a complete session
from a Python script.

## Before you begin

You need:

- LSQuant built and importable. See [Install LSQuant](install.md).
- The `mcp` Python package, in the same Python that imports LSQuant:

  ```bash
  pip install mcp
  ```

  Example 18 was tested with `mcp` 2.3.0.
- An MCP client. This guide uses Claude Code and Claude Desktop; other clients
  take the same command and environment.

## How the server works

The client starts the server as a child process and exchanges messages with it
only over standard input and output, so the server opens no network port.

- The server checks explicit run-file input paths, connector path parameters,
  generator scripts and output directories against one directory, the **root**.
  It resolves symbolic links and refuses paths outside the root before starting
  a run, with `Access denied: <key> … outside the server root`. Errors identify
  the supplied field without revealing a symbolic link's resolved target.
  The same check applies to run-file, inspection and reconstruction targets.
- The `run` tool returns at once with a `run_id`. The calculation runs in a
  separate `python -m lsquant run` process, which inherits the environment of
  the server, including `PYTHONPATH` and `OMP_NUM_THREADS`.
- Relative run-file inputs are tried against the root, then the run file's
  directory if absent at the root. DeepH's `hamiltonian_file` is relative to
  `model.dir`; GPAW's `dipole_cache` and `output.dir` are relative to the root.
  The default `output.dir` is `"."`, the root itself. This is allowed for
  compatibility; set `output.dir` for every run to keep outputs separate.

## Choose the root

The server uses the first of these settings that is present:

1. the option `--root DIR`;
2. the environment variable `LSQUANT_MCP_ROOT`;
3. if neither is set, the LSQuant source tree, two levels above `python/lsquant/`.

Use a working directory that holds your models and outputs, rather than your
home directory. For example, keep a project directory such as
`/path/to/projects/bi-shc` with the SIESTA files in `siesta/` and the runs in
`out/`, and pass it with `--root`. Put everything the assistant needs to read
inside the root, as described in [How the server works](#how-the-server-works).

## Test the server by hand

Run the server once in a terminal before you configure a client:

```bash
export PYTHONPATH="/path/to/lsquant/build/golden/lsquant_core:/path/to/lsquant/python"
python -m lsquant.mcp_server --root /path/to/projects/bi-shc
```

The server prints nothing and waits for MCP messages on standard input. Press
Ctrl+C to stop it. If it exits at once with `ModuleNotFoundError`, fix
`PYTHONPATH` or install `mcp`.

## Connect Claude Code

1. Register the server. Give the absolute path of the Python interpreter that
   imports LSQuant, and pass the environment with `-e`:

   ```bash
   claude mcp add lsquant \
       -e PYTHONPATH=/path/to/lsquant/build/golden/lsquant_core:/path/to/lsquant/python \
       -e OMP_NUM_THREADS=4 \
       -- /path/to/python -m lsquant.mcp_server --root /path/to/projects/bi-shc
   ```

   By default the server is registered for the current project only. Add
   `--scope user` to use it in every project, or `--scope project` to write it
   into a `.mcp.json` file that you can share with your group.

2. Check that Claude Code can start the server:

   ```bash
   claude mcp list
   ```

3. Start Claude Code and ask, for example, "Describe the LSQuant word
   `Tr[VX delta(E-H) VY dG-(E)/dE]`." The assistant calls the
   `words_describe` tool.

## Connect Claude Desktop

Claude Desktop reads its servers from `claude_desktop_config.json`. Open it
from the settings of the app (Developer, Edit Config), and add an entry under
`mcpServers`:

```json
{
  "mcpServers": {
    "lsquant": {
      "command": "/path/to/python",
      "args": ["-m", "lsquant.mcp_server", "--root", "/path/to/projects/bi-shc"],
      "env": {
        "PYTHONPATH": "/path/to/lsquant/build/golden/lsquant_core:/path/to/lsquant/python",
        "OMP_NUM_THREADS": "4"
      }
    }
  }
}
```

Restart Claude Desktop after you save the file.

Claude Desktop starts the server on the computer where it runs. If LSQuant is
installed on a remote machine, you can let the client start the server over
`ssh`. This works because the server communicates only through standard input
and output (see [How the server works](#how-the-server-works)). It needs a
passwordless, key-based login:

```json
{
  "mcpServers": {
    "lsquant": {
      "command": "ssh",
      "args": ["user@workstation",
               "PYTHONPATH=/path/to/lsquant/build/golden/lsquant_core:/path/to/lsquant/python OMP_NUM_THREADS=4 /path/to/python -m lsquant.mcp_server --root /path/to/projects/bi-shc"]
    }
  }
}
```

Run the server on a workstation or in an interactive session on a compute
node, not on a cluster login node. Every run that the server starts executes
on the same machine as the server, so a server on a login node runs the
calculations there.

## Set the environment

Put these variables in the `env` block of the client, or in the `-e` options
of `claude mcp add`:

| Variable | Why |
|---|---|
| `PYTHONPATH` | The server and every run it starts must import `lsquant` and `lsquant._core`. |
| `OMP_NUM_THREADS` | Number of threads of each run. If you leave it unset, `import lsquant` sets it to every physical core the process may use, and each run of the assistant then takes the whole machine. |
| `LSQUANT_MCP_ROOT` | The root, when `--root` is absent. |

Each `run` call starts a new process, and the server does not limit how many
run at the same time. If the assistant starts several runs, they share the
cores; set `OMP_NUM_THREADS` low enough for that.

## Work with the assistant

A typical session:

1. Ask the assistant to describe the word you want, to check that LSQuant
   accepts it (`words_describe`).
2. Give the assistant a run file in the root, or let it write one, and ask it
   to run it (`run`).
3. Ask for the status (`status`). The log is in
   `<output_dir>/run_<run_id>.log`.
4. Ask it to inspect the output directory (`inspect`) and to reconstruct
   with another broadening (`reconstruct`).

The server keeps the run identifiers in memory. After the client restarts the
server, `status` reports only the runs started since the restart. The output
directories and logs of earlier runs remain on disk, and `inspect` still
reads them.

When you request a formula for a directory containing both rank-1 (DOS) and
rank-2 (conductivity) moments, `reconstruct` skips incompatible files and
returns their paths and reasons in `skipped_files`. If no file admits the
formula, it reports an error before writing curves. Without an explicit
formula, each file uses its own default; see example 18.

The run-file path guard includes `model.file`, `model.hr_file`, string-valued
`model.positions` and `model.atom_map`, and the connector paths `overlap_file`,
`stru_file`, `orbital_dir`, `position_file`, `dipole_cache` (read and written),
and `hamiltonian_file` in `model.params`. Inline arrays are unaffected.
These checks validate explicit paths before dispatch; they do not sandbox
connector code or prevent concurrent filesystem changes. `status` accepts a
run identifier, not a filesystem path.

## See also

- [MCP tool reference](../reference/mcp-tools.md)
- [Example 18: an MCP session](../../examples/18_mcp_session/README.md)
- [How we run LSQuant](run-lsquant.md)
- [Run-file reference](../reference/run-file.md)
