# mlstdb

[![Tests](https://github.com/MDU-PHL/mlstdb/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/MDU-PHL/mlstdb/actions/workflows/test.yml)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/MDU-PHL/mlstdb)](https://github.com/MDU-PHL/mlstdb/releases)
[![PyPI - Version](https://img.shields.io/pypi/v/mlstdb.svg)](https://pypi.org/project/mlstdb)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mlstdb.svg)](https://pypi.org/project/mlstdb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/mlstdb/badges/version.svg)](https://anaconda.org/bioconda/mlstdb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/mlstdb/badges/license.svg)](https://anaconda.org/bioconda/mlstdb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/mlstdb/badges/downloads.svg)](https://anaconda.org/bioconda/mlstdb)

Keep your [`mlst`](https://github.com/tseemann/mlst) databases up to date. `mlstdb` handles OAuth authentication with [PubMLST](https://pubmlst.org/) and [BIGSdb Pasteur](https://bigsdb.pasteur.fr/) so you can download the latest MLST schemes and build a BLAST database, in two commands.

**[Full Documentation](https://MDU-PHL.github.io/mlstdb)**

## Install

```sh
conda create -n mlst -c bioconda mlst
conda activate mlst
pip install mlstdb
```

<details>
<summary>Other installation methods</summary>

```sh
# From bioconda (include conda-forge for dependencies)
conda install -c conda-forge -c bioconda mlstdb

# Or install both tools together
conda create -n mlst -c conda-forge -c bioconda mlst mlstdb

# From PyPI only
pip install mlstdb
```

</details>

## Quick Start

**1. Register with each database (one-time setup):**

```sh
mlstdb connect --db pubmlst
mlstdb connect --db pasteur
```

This opens a browser for OAuth registration. Follow the prompts to authorise `mlstdb`.

**2. Download schemes and build the BLAST database:**

```sh
mlstdb update
```

This downloads the curated MLST schemes from both PubMLST and Pasteur and creates a BLAST database.

**3. Use with `mlst`:**

```sh
mlst --blastdb blast/mlst.fa --datadir pubmlst your_assembly.fasta
```

That's it. For advanced scheme exploration, custom filtering, and detailed option reference, see the [full documentation](https://MDU-PHL.github.io/mlstdb).

## Caution

- **Back up** your existing MLST databases before running updates.
- **Curate** your scheme list before updating — not all schemes are validated for the `mlst` tool.
- The `mlst` tool is designed for **bacterial species only**.

## Acknowledgements

Built upon the work of:

- [BIGSdb_downloader](https://github.com/kjolley/BIGSdb_downloader) by Keith Jolley
- [pyMLST](https://github.com/bvalot/pyMLST) by Benoit Valot

## License

GPL v3. See [LICENSE.txt](LICENSE.txt) for details.
