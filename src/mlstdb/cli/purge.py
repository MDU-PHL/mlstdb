import sys
import click
from pathlib import Path

from mlstdb.core.download import create_blast_db
from mlstdb.core.purge import (
    parse_purge_config,
    purge_scheme,
    purge_st,
    purge_allele,
    purge_st_and_allele,
)
from mlstdb.utils import error, info, success


def _run_single_purge(scheme, st, allele, force, verbose, directory):
    """Execute a single purge operation (no BLAST rebuild).

    Returns True if a change was made, False otherwise.
    """
    scheme_dir = Path(directory) / scheme

    if not scheme_dir.is_dir():
        error(f"Scheme directory not found: {scheme_dir}")
        return False

    # Parse allele string if provided
    locus = None
    allele_num = None
    if allele is not None:
        if ":" not in allele:
            error("Invalid allele format. Expected locus:number (e.g. aroC:3).")
            return False
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
                return False
        purge_scheme(str(scheme_dir), verbose=verbose)
        return True

    elif st is not None and allele is None:
        # Mode 2: Purge a single ST
        return purge_st(
            str(scheme_dir), scheme, st,
            force=force, verbose=verbose,
        )

    elif st is None and allele is not None:
        # Mode 3: Purge an allele (and affected STs)
        return purge_allele(
            str(scheme_dir), scheme, locus, allele_num,
            force=force, verbose=verbose,
        )

    else:
        # Mode 4: Purge a specific ST and check a specific allele
        return purge_st_and_allele(
            str(scheme_dir), scheme, st, locus, allele_num,
            force=force, verbose=verbose,
        )


def _run_config_purge(config_path, force, verbose, directory, blast_directory):
    """Execute batch purge from a YAML config file.

    CLI flags override global options in the config.
    """
    try:
        entries, global_opts = parse_purge_config(config_path)
    except (ValueError, Exception) as e:
        error(f"Failed to parse config file: {e}")
        sys.exit(1)

    # CLI flags take precedence; fall back to config globals, then defaults
    eff_force = force or global_opts["force"]
    eff_verbose = verbose or global_opts["verbose"]
    eff_directory = directory if directory != "pubmlst" else (global_opts["directory"] or "pubmlst")
    eff_blast_dir = blast_directory or global_opts["blast_directory"]

    total_ops = sum(
        max(len(e["st"]) + len(e["alleles"]), 1) for e in entries
    )
    info(f"Config loaded: {len(entries)} scheme(s), {total_ops} operation(s).")

    changes = 0
    errors = 0

    for entry in entries:
        scheme = entry["scheme"]
        sts = entry["st"]
        alleles = entry["alleles"]

        if not sts and not alleles:
            # Purge entire scheme
            info(f"\n--- Purging scheme: {scheme} ---")
            ok = _run_single_purge(scheme, None, None, eff_force, eff_verbose, eff_directory)
            if ok:
                changes += 1
            else:
                errors += 1
            continue

        # Process STs first
        for st in sts:
            info(f"\n--- Purging ST {st} from {scheme} ---")
            ok = _run_single_purge(scheme, st, None, eff_force, eff_verbose, eff_directory)
            if ok:
                changes += 1
            else:
                errors += 1

        # Then process alleles
        for allele in alleles:
            info(f"\n--- Purging allele {allele} from {scheme} ---")
            ok = _run_single_purge(scheme, None, allele, eff_force, eff_verbose, eff_directory)
            if ok:
                changes += 1
            else:
                errors += 1

    info(f"\nBatch summary: {changes} successful, {errors} failed.")

    if changes == 0:
        error("No changes were made. BLAST database rebuild skipped.")
        sys.exit(1)

    # Single BLAST rebuild at the end
    info("\nRebuilding BLAST database...")
    create_blast_db(eff_directory, eff_blast_dir, eff_verbose)
    success("Batch purge completed successfully!")


@click.command()
@click.help_option("-h", "--help")
@click.option(
    "--scheme", "-s", default=None,
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
    "--config", "-c", "config_path", default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to a YAML config file for batch purging across multiple schemes.",
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
def purge(scheme: str, st: str, allele: str, config_path: str, force: bool,
          verbose: bool, directory: str, blast_directory: str):
    """
    Purge MLST scheme data and rebuild BLAST database.

    Remove an entire scheme, a specific ST, or a specific allele from the local
    database. The BLAST database is rebuilt afterwards.

    Use --config to batch-purge across multiple schemes from a YAML file.
    The BLAST database is rebuilt only once at the end.

    \b
    Examples:
      mlstdb purge -s salmonella              # Remove entire scheme
      mlstdb purge -s salmonella --st 3       # Remove ST 3
      mlstdb purge -s salmonella -a aroC:1    # Remove allele aroC_1
      mlstdb purge -s salmonella --st 3 -f    # Force, no prompts
      mlstdb purge -c purge_config.yaml       # Batch purge from config
    """
    try:
        # Validate option combinations
        if config_path and scheme:
            error("Cannot use --config together with --scheme. Use one or the other.")
            sys.exit(1)
        if not config_path and not scheme:
            error("Either --scheme or --config is required.")
            sys.exit(1)
        if config_path and (st or allele):
            error("Cannot use --config together with --st or --allele. "
                  "Define these in the config file instead.")
            sys.exit(1)

        # Batch mode
        if config_path:
            _run_config_purge(config_path, force, verbose, directory, blast_directory)
            return

        # Single mode
        changed = _run_single_purge(scheme, st, allele, force, verbose, directory)
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
