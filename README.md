# mlstdb

[![PyPI - Version](https://img.shields.io/pypi/v/mlstdb.svg)](https://pypi.org/project/mlstdb)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mlstdb.svg)](https://pypi.org/project/mlstdb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/mlstdb/badges/version.svg)](https://anaconda.org/bioconda/mlstdb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/mlstdb/badges/license.svg)](https://anaconda.org/bioconda/mlstdb)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/mlstdb/badges/downloads.svg)](https://anaconda.org/bioconda/mlstdb)

`mlstdb` is a Python package to update and manage the MLST database for the `mlst` tool using the PubMLST and BIGSdb Pasteur APIs. It is written to handle the OAuth2 authentication process that's required to access up-to-date MLST schemes available on these databases. This tool allows user to fetch MLST schemes, filter the schemes, and update the MLST database for the `mlst` tool.

-----

## Table of Contents

- [mlstdb](#mlstdb)
  - [Table of Contents](#table-of-contents)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Final Steps](#final-steps)
  - [License](#license)

## Prerequisites

Should install `mlst` for the use of this tool.  

## Installation

From bioconda: 
```sh
conda install -c bioconda mlstdb
```

From PyPI:
```sh
pip install mlstdb
```

Recommended way to install is from bioconda with `mlst` tool. 

```sh
conda create -n mlst -c bioconda mlst mlstdb
```

## Usage

`mlstdb` uses a simple two step process to update the MLST database for the `mlst` tool. It has two main subcommands: `fetch` and `update`.

1. **Fetch MLST schemes**

```sh
mlstdb fetch --help
```

```console
Usage: mlstdb fetch [OPTIONS]

  BIGSdb Scheme Fetcher Tool

  This tool downloads MLST scheme information from BIGSdb databases. It will
  automatically handle authentication and save the results.

Options:
  -h, --help                  Show this message and exit.
  -d, --db [pubmlst|pasteur]  Database to use (pubmlst or pasteur)
  -e, --exclude TEXT          Scheme name must not include provided term
                              (default: cgMLST)
  -m, --match TEXT            Scheme name must include provided term (default:
                              MLST)
  -s, --scheme-uris TEXT      Optional: Path to custom scheme_uris.tab file
  -f, --filter TEXT           Filter species or schemes using a wildcard
                              pattern
  -r, --resume                Resume processing from where it stopped
  -v, --verbose               Enable verbose logging for debugging
```

Use the `fetch` command to download MLST schemes from the BIGSdb databases. The `--db` argument specifies the database to use, which can be either `pubmlst` or `pasteur`. The `--exclude` and `--match` arguments can be used to filter the schemes based on the scheme name. The `--scheme-uris` argument can be used to provide a custom scheme URIs file. The `--filter` argument can be used to filter species or schemes using a wildcard pattern. The `--resume` flag can be used to resume processing from where it stopped. The `--verbose` flag can be used to enable verbose logging for debugging. This will create a `mlst_schemes_<db>.txt` file with the MLST schemes.

We can just use `mlstdb fetch` to download the MLST schemes from the BIGSdb databases. The command will prompt for the `db` (either `pubmlst` or `pasteur`) to fetch. If the registration is not done, it will prompt the user to register the client credentials. This will save the client credentials to the `~/.config/mlstdb` directory.

In cases where the tool does not find an appropriate scheme name, it will prompt the user to either set the missing schemes as 'missing' or auto-generate them. The user can choose the appropriate option as they are prompted.

<details>
<summary>Auto extraction of scheme?🤔</summary>

First, the script automatically tries to extract the scheme names from the `dbases.sh` file. If the scheme name is not found, it will prompt the user to either print `missing` in the output file or automatically create a scheme name based on the URL. For eg, for URL `https://rest.pubmlst.org/db/pubmlst_borrelia_seqdef/schemes/1`, the scheme name will be `borrelia`. If there are multiple schemes, it will append a number to the scheme name. For eg, for URLs `https://rest.pubmlst.org/db/pubmlst_chlamydiales_seqdef/schemes/38` and `https://rest.pubmlst.org/db/pubmlst_chlamydiales_seqdef/schemes/41`, the scheme names will be `chlamydiales_38` and `chlamydiales_41` respectively.

</details>


The script offers feature to filter for particular species/schemes. It is recommended to run with filter option and thus, download only the required schemes so as not to tamper with the existing DBs and schemes.

**📝Important**: `mlst` tool is designed for typing bacterial species only. Please make sure to filter the non-bacterial schemes from your schemes file.


2. **Update MLST database**

```sh
mlstdb update --help
```

```console
Usage: mlstdb update [OPTIONS]

  Update MLST schemes and create BLAST database.

  Downloads MLST schemes from the specified input file and creates a BLAST
  database from the downloaded sequences. Authentication tokens should be set
  up using fetch.py.

Options:
  -h, --help                  Show this message and exit.
  -i, --input TEXT            Path to mlst_schemes_<db>.tab containing MLST
                              scheme URLs  [required]
  -d, --directory TEXT        Directory to save the downloaded MLST schemes
                              (default: pubmlst)
  -b, --blast-directory TEXT  Directory for BLAST database (default: blast)
  -v, --verbose               Enable verbose logging for debugging
```

Use the `update` command to update the MLST database and create a BLAST database. The `--input` argument specifies the path to the `mlst_schemes_<db>.tab` file containing MLST scheme URLs. The `--directory` argument specifies the directory to save the downloaded MLST schemes. The `--blast-directory` argument specifies the directory for the BLAST database. The `--verbose` flag can be used to enable verbose logging for debugging.

We can prepare a custom `mlst_schemes_<db>.tab` file with headers `database	species	scheme_description	scheme	URI`
and use `mlstdb update` to update the MLST database for select species and schemes. This will automatically create a BLAST database from the downloaded sequences.

## Final Steps
After running all scripts, verify the database setup by running the `mlst` tool with the updated database:
```bash
mlst --blastdb <path_to_blast/mlst.fa> --datadir <path_to_pubmlst_dir>
```


## License

`mlstdb` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.


For additional support, please raise an issue.