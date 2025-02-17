# SPDX-FileCopyrightText: 2025-present Himal Shrestha <stha.himal2007@gmail.com>
#
# SPDX-License-Identifier: MIT

import click
from mlstdb.cli.fetch import fetch
from mlstdb.cli.update import update
from mlstdb.__about__ import __version__

@click.group()
@click.version_option(version=__version__)

@click.help_option('-h', '--help')
def mlstdb():
    """MLST Database Management Tool"""
    pass

mlstdb.add_command(fetch)
mlstdb.add_command(update)