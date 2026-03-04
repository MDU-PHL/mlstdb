# Installation

## Prerequisites

You need the `mlst` tool installed before using `mlstdb`.

## Recommended (pip + conda)

```sh
conda create -n mlst -c bioconda mlst
conda activate mlst
pip install mlstdb
```

## Via Bioconda

```sh
conda install -c conda-forge -c bioconda mlstdb
```

> **Note:** If you see `nothing provides rauth >=0.7.3`, add `-c conda-forge`.

## Via PyPI only

```sh
pip install mlstdb
```
