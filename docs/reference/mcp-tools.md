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
| `--root DIR` | Root directory for explicit input and output paths; see [Root directory](#root-directory). |

Requirements:

- The `mcp` package: `pip install mcp`. It is needed only for the server.
- `PYTHONPATH` must contain `build/<preset>/lsquant_core` and `python/`, as
  for the command-line interface.

## Root directory

The server checks run-file paths, all explicit model input paths, connector
path parameters, generator scripts and output directories after resolving
symbolic links. Escapes are refused before dispatch, with the offending key
and supplied value in the error, without revealing resolved outside paths.
Inspection and reconstruction targets receive the same check; `status` takes
only a run identifier. See [Path resolution](run-file.md#path-resolution) for
input fallback and connector-specific bases. The default `output.dir: "."`
is allowed and writes at the root. These checks do not sandbox connector code
or protect against concurrent filesystem changes.

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
