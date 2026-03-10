import sys
import click
import importlib.resources
from pathlib import Path

from mlstdb.core.auth import get_client_credentials, retrieve_session_token
from mlstdb.core.download import get_mlst_files, create_blast_db
from mlstdb.core.config import check_dir
from mlstdb.utils import error, success, info
from tqdm import tqdm


@click.command()
@click.help_option('-h', '--help')
@click.option('--input', '-i', default=None,
              help='Optional: Path to mlst_schemes_<db>.tab containing MLST scheme URLs')
@click.option('--directory', '-d', default='pubmlst',
              help='Directory to save the downloaded MLST schemes (default: pubmlst)')
@click.option('--blast-directory', '-b',
              help='Directory for BLAST database (default: blast)')
@click.option('--verbose', '-v', is_flag=True,
              help='Enable verbose logging for debugging')
@click.option('--no-auth', 'no_auth', is_flag=True, default=False,
              help='Skip OAuth authentication and use unauthenticated requests. Schemes requiring authentication will be skipped.')
def update(input: str, directory: str, blast_directory: str, verbose: bool, no_auth: bool):
    """
    Update MLST schemes and create BLAST database.

    Downloads MLST schemes from PubMLST and Pasteur databases using the curated scheme list included in the package or using a user-specified input file and creates a BLAST database from the downloaded sequences.
    
    ---
    
    Please ensure you have:
    
    1. Successful connection to both PubMLST and Pasteur databases via 'mlstdb connect'
    
    2. Registered to all the databases you wish to update.
    
    ---
    """
    try:
        if input is None:
            with importlib.resources.path('mlstdb.data', 'mlst_schemes_all.tab') as default_path:
                input = str(default_path)

        with open(input, 'r') as f:
            next(f)  # Skip header
            lines = f.readlines()

        check_dir(directory)

        skipped_schemes = []
        download_success = False

        for line in tqdm(lines, desc="Downloading MLST schemes", unit="scheme"):
            parts = line.strip().split('\t')
            if len(parts) != 5:
                error(f"Skipping invalid line: {line}")
                continue

            database, species, scheme_desc, scheme, url = parts

            if no_auth:
                client_key = client_secret = session_token = session_secret = None
            else:
                try:
                    client_key, client_secret = get_client_credentials(database.lower())
                    session_token, session_secret = retrieve_session_token(database.lower())
                except ValueError as ve:
                    error(f"Error with credentials for {database}: {ve}")
                    skipped_schemes.append((scheme, database))
                    continue

                if not session_token or not session_secret:
                    error(f"No valid session token found for {database}.")
                    skipped_schemes.append((scheme, database))
                    continue

            scheme_dir = Path(directory) / scheme
            check_dir(str(scheme_dir))

            try:
                get_mlst_files(url, str(scheme_dir), client_key, client_secret,
                               session_token, session_secret, scheme, verbose=verbose,
                               no_auth=no_auth)
                success(f"Successfully downloaded scheme: {scheme}")
                download_success = True
            except Exception as e:
                if '401 Client Error' in str(e) or '403 Client Error' in str(e):
                    if no_auth:
                        error(f"Scheme [{scheme}] requires authentication and cannot be downloaded with --no-auth.")
                        info(f"Run without --no-auth to download authenticated schemes.")
                    else:
                        error(f"Authentication error for [{scheme}]: {e}")
                        error(f"Please make sure you have registered to [{scheme}] scheme within the {database} database.")
                else:
                    error(f"Error downloading scheme [{scheme}]: {e}")
                    if not no_auth:
                        info(f"Please make sure you have registered to [{scheme}] scheme within the {database} database.")
                skipped_schemes.append((scheme, database))

        if skipped_schemes:
            click.secho(
                "\nThe following schemes were not downloaded because you are not registered to them:",
                fg="yellow"
            )
            for scheme_name, db_name in skipped_schemes:
                click.secho(f"  - {scheme_name} from {db_name}", fg="yellow")
            click.secho(
                "\nFor the complete update, please ensure you have registered to the above schemes "
                "and then run this command again.", fg="yellow"
            )
            click.secho(
                "If you choose to continue, the BLAST database will be created with only the schemes "
                "you have downloaded.", fg="yellow"
            )

            if not download_success:
                error("\nNo schemes were successfully downloaded. BLAST database creation skipped.")
                sys.exit(1)

            if not click.confirm(
                "\nAre you sure you want to continue with BLAST database creation "
                "with the schemes you are registered to and have downloaded? (y/n)",
                default=False,
                prompt_suffix=" "
            ):
                info("BLAST database creation cancelled.")
                sys.exit(0)

        if not [d for d in Path(directory).iterdir() if d.is_dir()]:
            error("\nNo schemes were successfully downloaded. BLAST database creation skipped.")
            sys.exit(1)

        info("\nCreating BLAST database from downloaded MLST schemes...")
        create_blast_db(directory, blast_directory, verbose)
        success("Update completed successfully!")

    except FileNotFoundError:
        error(f"Input file not found: {input}")
        info("Please run 'mlstdb fetch' first to generate the scheme list file.")
        sys.exit(1)
    except Exception as e:
        error(f"An error occurred: {e}")
        if verbose:
            import traceback
            error(traceback.format_exc())
        sys.exit(1)