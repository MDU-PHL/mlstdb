# mlstdb

!!!!!!IN CONSTRUCTION!!!!!

`mlstdb` is a Python package to update and manage the MLST database for the `mlst` tool using the PubMLST and BIGSdb Pasteur APIs.

[![Tests](https://github.com/MDU-PHL/mlstdb/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/MDU-PHL/mlstdb/actions/workflows/test.yml)
[![PyPI](https://img.shields.io/pypi/v/mlstdb.svg)](https://pypi.org/project/mlstdb)
[![Bioconda](https://anaconda.org/bioconda/mlstdb/badges/version.svg)](https://anaconda.org/bioconda/mlstdb)

## Overview

mlstdb handles OAuth2 authentication to fetch up-to-date MLST schemes from PubMLST and BIGSdb Pasteur, and updates your local `mlst` tool database.