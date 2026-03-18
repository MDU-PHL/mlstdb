# SPDX-FileCopyrightText: 2025-present Himal Shrestha <stha.himal2007@gmail.com>
#
# SPDX-License-Identifier: GPL-3.0-or-later

import click

from mlstdb.__about__ import __version__
from mlstdb.cli.connect import connect
from mlstdb.cli.fetch import fetch
from mlstdb.cli.purge import purge
from mlstdb.cli.update import update

@click.group(add_help_option=False)
@click.version_option(version=__version__)
@click.help_option("-h", "--help")
def mlstdb():
    """MLST Database Management Tool.

    \b
    Recommended workflow:
    1. mlstdb connect --db pubmlst/pasteur    # Register credentials
    2. mlstdb update --db pubmlst             # Update database
    3. mlstdb fetch                           # [ADVANCED] Explore custom schemes
    
    Most users only need to use 'connect' and 'update' commands.
    
    """
    pass


mlstdb.add_command(connect)
mlstdb.add_command(update)
mlstdb.add_command(fetch)
mlstdb.add_command(purge)