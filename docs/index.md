# mlstdb

[![Tests](https://github.com/MDU-PHL/mlstdb/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/MDU-PHL/mlstdb/actions/workflows/test.yml)
[![PyPI](https://img.shields.io/pypi/v/mlstdb.svg)](https://pypi.org/project/mlstdb)
[![Bioconda](https://anaconda.org/bioconda/mlstdb/badges/version.svg)](https://anaconda.org/bioconda/mlstdb)

**Keep your [`mlst`](https://github.com/tseemann/mlst) databases up to date with.. just two commands**

`mlstdb` handles the OAuth authentication required to access [PubMLST](https://pubmlst.org/) and [BIGSdb Pasteur](https://bigsdb.pasteur.fr/) APIs, downloads the latest MLST schemes, and builds a ready-to-use BLAST database for the `mlst` tool.

---

## Why `mlstdb`?

The `mlst` tool comes with a bundled database, but MLST schemes are continuously updated on PubMLST and Pasteur. Keeping your local database current requires authentication setup (OAuth2) and downloading files. `mlstdb` handles all of that with additional features.

**What it does:**

- Handles OAuth registration and token management for PubMLST and Pasteur
- Downloads allele sequences and ST profiles for curated MLST schemes
- Builds the BLAST database that `mlst` needs
- Supports parallel downloads, resume on failure, and custom scheme lists

---

## Quick Start

```sh
# 1. Install
conda create -n mlst -c bioconda mlst && conda activate mlst
pip install mlstdb

# 2. Connect to databases (one-time setup)
mlstdb connect --db pubmlst
mlstdb connect --db pasteur

# 3. Download schemes and build BLAST database
mlstdb update

# 4. Use with mlst
mlst --blastdb blast/mlst.fa --datadir pubmlst your_assembly.fasta
```

See the [Getting Started](getting-started.md) guide for a detailed walkthrough.

---

## Commands at a Glance

| Command | Purpose |
|---------|---------|
| `mlstdb connect` | Register OAuth credentials with PubMLST or Pasteur |
| `mlstdb update` | Download schemes and build the BLAST database |
| `mlstdb fetch` | *(Advanced)* Explore and filter all available schemes |

Most users only need `connect` and `update`. See the [Usage Overview](usage/overview.md) for details.

---

## Next Steps

- [Installation](installation.md) — All installation methods
- [Getting Started](getting-started.md) — End-to-end tutorial
- [Usage Overview](usage/overview.md) — How the commands work together
