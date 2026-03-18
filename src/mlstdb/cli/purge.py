import sys
import click
from pathlib import Path

from mlstdb.core.download import create_blast_db
from mlstdb.core.purge import (
    purge_scheme,
    purge_st,
    purge_allele,
    purge_st_and_allele,
)
from mlstdb.utils import error, info, success


@click.command()
@click.help_option("-h", "--help")
@click.option(
    "--scheme", "-s", required=True,
    help="Scheme name to purge (e.g. salmonella). Must match a directory under the data directory.",
)
@click.option(
    "--st", default=None,
    help="ST number to purge from the scheme.",
)
@click.option(
    "--allele", "-a", default=None,
    help="Allele to purge, format locus:number (e.g. aroC:3).",
)
@click.option(
    "--force", "-f", is_flag=True,
    help="Skip confirmation prompts.",
)
@click.option(
    "--verbose", "-v", is_flag=True,
    help="Enable verbose logging for debugging.",
)
@click.option(
    "--directory", "-d", default="pubmlst",
    help="Directory containing the downloaded MLST schemes (default: pubmlst).",
)
@click.option(
    "--blast-directory", "-b", default=None,
    help="Directory for BLAST database (default: blast).",
)
def purge(scheme: str, st: str, allele: str, force: bool, verbose: bool,
          directory: str, blast_directory: str):
    """
    Purge MLST scheme data.

    Remove an entire scheme, a specific ST, or a specific allele from the local
    database. The BLAST database is rebuilt afterwards.

    \b
    Examples:
      mlstdb purge -s salmonella              # Remove entire scheme
      mlstdb purge -s salmonella --st 3       # Remove ST 3
      mlstdb purge -s salmonella -a aroC:1    # Remove allele aroC_1
      mlstdb purge -s salmonella --st 3 -f    # Force, no prompts
    """
    try:
        scheme_dir = Path(directory) / scheme

        if not scheme_dir.is_dir():
            error(f"Scheme directory not found: {scheme_dir}")
            sys.exit(1)

        # Parse --allele if provided
        locus = None
        allele_num = None
        if allele is not None:
            if ":" not in allele:
                error("Invalid allele format. Expected locus:number (e.g. aroC:3).")
                sys.exit(1)
            locus, allele_num = allele.split(":", 1)

        # Determine purge mode
        if st is None and allele is None:
            # Mode 1: Purge entire scheme
            if not force:
                if not click.confirm(
                    f"This will delete the entire scheme directory '{scheme_dir}'. Continue?",
                    default=False,
                ):
                    info("Purge cancelled.")
                    sys.exit(0)
            purge_scheme(str(scheme_dir), verbose=verbose)

        elif st is not None and allele is None:
            # Mode 2: Purge a single ST
            changed = purge_st(
                str(scheme_dir), scheme, st,
                force=force, verbose=verbose,
            )
            if not changed:
                sys.exit(1)

        elif st is None and allele is not None:
            # Mode 3: Purge an allele (and affected STs)
            changed = purge_allele(
                str(scheme_dir), scheme, locus, allele_num,
                force=force, verbose=verbose,
            )
            if not changed:
                sys.exit(1)

        else:
            # Mode 4: Purge a specific ST and check a specific allele
            changed = purge_st_and_allele(
                str(scheme_dir), scheme, st, locus, allele_num,
                force=force, verbose=verbose,
            )
            if not changed:
                sys.exit(1)

        # Rebuild BLAST database
        info("\nRebuilding BLAST database...")
        create_blast_db(directory, blast_directory, verbose)
        success("Purge completed successfully!")

    except Exception as e:
        error(f"An error occurred: {e}")
        if verbose:
            import traceback
            error(traceback.format_exc())
        sys.exit(1)
