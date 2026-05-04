# Installation

## Prerequisites

`mlstdb` is designed to work alongside the [`mlst`](https://github.com/tseemann/mlst) tool. You'll need:

- **Python** >= 3.8
- **mlst** — install via [bioconda](https://bioconda.github.io/)
- **BLAST+** — comes with `mlst`, or install separately via bioconda/apt

## Recommended Method

Create a conda environment with `mlst`, then install `mlstdb` via pip:

```sh
conda create -n mlst -c bioconda mlst
conda activate mlst
pip install mlstdb
```

This is the most reliable method and avoids dependency conflicts.

## Alternative Methods

=== "Bioconda"

    ```sh
    conda install -c conda-forge -c bioconda mlstdb
    ```

    !!! note
        Include the `-c conda-forge` channel to resolve all dependencies. If you see errors like `nothing provides rauth >=0.7.3`, this is the fix.

=== "Both tools together"

    ```sh
    conda create -n mlst -c conda-forge -c bioconda mlst mlstdb
    ```

=== "PyPI only"

    ```sh
    pip install mlstdb
    ```

    !!! warning
        You'll need to install `mlst` and BLAST+ separately if using pip only.

## Verify Installation

```sh
mlstdb --version
mlstdb --help
```

You should see the version number and the available commands (`connect`, `update`, `fetch`).

## Next Steps

Head to the [Getting Started](getting-started.md) guide to connect to the databases and download your first MLST schemes.
