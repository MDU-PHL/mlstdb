# mlstdb

[![Tests](https://github.com/MDU-PHL/mlstdb/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/MDU-PHL/mlstdb/actions/workflows/test.yml)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/MDU-PHL/mlstdb)](https://github.com/MDU-PHL/mlstdb/releases)
[![PyPI - Version](https://img.shields.io/pypi/v/mlstdb.svg)](https://pypi.org/project/mlstdb)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mlstdb.svg)](https://pypi.org/project/mlstdb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/mlstdb/badges/version.svg)](https://anaconda.org/bioconda/mlstdb)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
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

## Removing contaminated STs or alleles

Discovered a dodgy sequence type or allele in your local database? You can remove it without re-downloading the whole scheme:

```sh
# Remove a single ST (orphaned alleles are cleaned up automatically)
mlstdb purge --scheme salmonella --st 3

# Remove a specific allele (also removes any STs that reference it)
mlstdb purge --scheme salmonella --allele aroC:1
```

The BLAST database is rebuilt automatically after each purge. See the [purge documentation](https://MDU-PHL.github.io/mlstdb/usage/purge/) for the full reference.

## Caution

- **Back up** your existing MLST databases before running updates.
- If using `mlstdb fetch` to build a custom scheme list, **double-check** that all schemes are compatible with the `mlst` tool. Not all schemes are validated for use with `mlst`. The `mlst` tool is designed for **bacterial species only**.
- `mlstdb purge` permanently modifies your local database. Take a backup of your `pubmlst/` directory before purging, especially when using `--force`.

## Acknowledgements

Built upon the work of:

- [BIGSdb_downloader](https://github.com/kjolley/BIGSdb_downloader) by Keith Jolley
- [pyMLST](https://github.com/bvalot/pyMLST) by Benoit Valot

## License

`mlstdb` was previously licensed under MIT. As of version 0.1.7, it is licensed under GPL v3. Original MIT‑licensed code is preserved and attributed according to MIT terms.

For additional support, please raise an issue.
