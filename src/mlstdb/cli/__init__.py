# SPDX-FileCopyrightText: 2025-present Himal Shrestha <stha.himal2007@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

import click
from mlstdb.cli.fetch import fetch
from mlstdb.cli.update import update
from mlstdb.cli.inspect import inspect
from mlstdb.cli.compare import compare
from mlstdb.__about__ import __version__

@click.group()
@click.version_option(version=__version__)

@click.help_option('-h', '--help')
def mlstdb():
    """MLST Database Management Tool"""
    pass

mlstdb.add_command(fetch)
mlstdb.add_command(update)
mlstdb.add_command(inspect)
mlstdb.add_command(compare)