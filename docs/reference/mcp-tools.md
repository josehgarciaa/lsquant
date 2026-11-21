# MCP server

`python -m lsquant.mcp_server` runs a [Model Context Protocol](https://modelcontextprotocol.io)
(MCP) server over standard input and output. An MCP client, such as an AI
assistant, can then parse words, start runs and read results through tool
calls. The server is built on `FastMCP` from the `mcp` Python package.

For a step-by-step setup, see [Connect an AI assistant](../how-to/mcp-server.md).

## Start the server

```text
python -m lsquant.mcp_server [--root DIR]
```

| Option | Meaning |
|---|---|
| `--root DIR` | Root directory. The server reads and writes only inside it, except for `model.hr_file`, `model.positions` and `model.atom_map`, which it does not check; see [Root directory](#root-directory). |

Requirements:

- The `mcp` package: `pip install mcp`. It is needed only for the server.
- `PYTHONPATH` must contain `build/<preset>/lsquant_core` and `python/`, as
  for the command-line interface.

## Root directory

The server checks run files, `model.dir`, `model.csr_dir`, generator scripts
and output directories. It refuses any of them that resolves outside the root
directory. It resolves symbolic links before the check, and it takes a relative
path relative to the root. It does not check `model.hr_file`,
`model.positions` or `model.atom_map`, so a run file can read those files from
outside the root (`python/lsquant/mcp_server.py:108–118`). See the
[known limitation](../how-to/mcp-server.md#work-with-the-assistant) in
Connect an AI assistant.

The root is chosen in this order:

1. `--root DIR`;
2. the environment variable `LSQUANT_MCP_ROOT`;
3. the repository root (two levels above `python/lsquant/`).

## Tools

All tools return a JSON string. If a tool fails, the client receives an error
with the message of the underlying exception.

### `words_describe(word)`

| Parameter | Type | Meaning |
|---|---|---|
| `word` | string | Observable word, for example `Tr[delta(E-H)]`. |

Returns the same JSON as `python -m lsquant words describe "<word>"`.

### `run(runfile_json)`

| Parameter | Type | Meaning |
|---|---|---|
| `runfile_json` | string | Either the run file itself as a JSON string (it starts with `{`), or the path of a run file inside the root. |

The tool validates the run file and checks every path against the root. It then
starts `python -m lsquant run` in a background process and returns at once:

```json
{"run_id": "3f9c0a1b2d4e", "log": "<output_dir>/run_3f9c0a1b2d4e.log", "output_dir": "<output_dir>"}
```

It writes the resolved run file to `<output_dir>/.mcp_run_<run_id>.json` and
the log of the run to `<output_dir>/run_<run_id>.log`.

### `status(run_id)`

| Parameter | Type | Meaning |
|---|---|---|
| `run_id` | string | The identifier that `run` returned. |

Returns `run_id`, `status` (`running`, `succeeded` or `failed`), `returncode`
(`null` while running), `log` and `output_dir`.

The server keeps run identifiers in memory, so after a restart `status` knows
only the runs started since then. The output directories and logs of earlier
runs stay on disk, and you can read them with `inspect`.

### `reconstruct(output_dir, formula=None, broadening_meV=None, fermi_energies=None)`

| Parameter | Type | Meaning |
|---|---|---|
| `output_dir` | string | A run output directory or a moment file inside the root. |
| `formula` | string, optional | Reconstruction formula, as in [`reconstruct`](cli.md#reconstruct). |
| `broadening_meV` | number, optional | Target broadening in meV. |
| `fermi_energies` | list of numbers, optional | Fermi energies in eV, for time correlations. |

Returns the same summary as `python -m lsquant reconstruct --json`. The curves
are written next to the moments.

### `inspect(path)`

| Parameter | Type | Meaning |
|---|---|---|
| `path` | string | A moment file, run file or output directory inside the root. |

Returns the same JSON as `python -m lsquant inspect <path>`.

## Client configuration

Most MCP clients read a JSON block of this form. Replace the paths with yours:

```json
{
  "mcpServers": {
    "lsquant": {
      "command": "python3",
      "args": ["-m", "lsquant.mcp_server"],
      "env": {
        "PYTHONPATH": "/path/to/lsquant/build/golden/lsquant_core:/path/to/lsquant/python",
        "LSQUANT_MCP_ROOT": "/path/to/lsquant",
        "OMP_NUM_THREADS": "16"
      }
    }
  }
}
```

`OMP_NUM_THREADS` is passed on to the background runs, because each run inherits
the environment of the server.
