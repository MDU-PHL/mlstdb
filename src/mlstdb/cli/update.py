import sys
import click
import importlib.resources
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from mlstdb.core.auth import get_client_credentials, retrieve_session_token
from mlstdb.core.download import get_mlst_files, create_blast_db, create_session
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
@click.option('--resume', '-r', is_flag=True,
              help='Resume from where it stopped, skipping already downloaded schemes')
@click.option('--threads', '-t', default=1, show_default=True,
              help='Number of parallel threads for downloading schemes, max: 4 to avoid too many simultaneous requests')
def update(input: str, directory: str, blast_directory: str, verbose: bool, no_auth: bool,
           resume: bool, threads: int):
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

        if no_auth:
            click.secho(
                "\nWarning: Using unauthenticated access. The database is only available up to "
                "2024-12-31. Schemes created after this date may not be available.",
                fg="yellow"
            )

        # Parse all scheme entries
        scheme_entries = []
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) != 5:
                error(f"Skipping invalid line: {line}")
                continue
            scheme_entries.append(parts)

        # Filter already-downloaded schemes when resuming
        if resume:
            original_count = len(scheme_entries)
            remaining = []
            for entry in scheme_entries:
                database, species, scheme_desc, scheme, url = entry
                scheme_dir = Path(directory) / scheme
                profiles_file = scheme_dir / f"{scheme}.txt"
                if profiles_file.exists():
                    if verbose:
                        info(f"Skipping already downloaded scheme: {scheme}")
                else:
                    remaining.append(entry)
            skipped = original_count - len(remaining)
            if skipped > 0:
                info(f"Resuming: skipping {skipped} already downloaded scheme(s), "
                     f"{len(remaining)} remaining")
            scheme_entries = remaining

        # Build reusable sessions per database type
        sessions = {}  # db_type -> session
        credentials = {}  # db_type -> (client_key, client_secret, session_token, session_secret)

        if no_auth:
            sessions['no_auth'] = create_session(no_auth=True)
        else:
            db_types_needed = set(entry[0].lower() for entry in scheme_entries)
            for db_type in db_types_needed:
                try:
                    client_key, client_secret = get_client_credentials(db_type)
                    session_token, session_secret = retrieve_session_token(db_type)
                except ValueError as ve:
                    error(f"Error with credentials for {db_type}: {ve}")
                    continue
                if not session_token or not session_secret:
                    error(f"No valid session token found for {db_type}.")
                    continue
                credentials[db_type] = (client_key, client_secret, session_token, session_secret)
                sessions[db_type] = create_session(client_key, client_secret,
                                                    session_token, session_secret)

        skipped_schemes = []
        download_success = False
        parallel = threads > 1

        def _download_scheme(entry):
            """Worker function for downloading a single scheme."""
            database, species, scheme_desc, scheme, url = entry
            db_key = 'no_auth' if no_auth else database.lower()

            if db_key not in sessions:
                return ('skip_cred', scheme, database)

            session = sessions[db_key]
            creds = credentials.get(database.lower(), (None, None, None, None))
            client_key, client_secret, session_token, session_secret = creds

            scheme_dir = Path(directory) / scheme
            check_dir(str(scheme_dir))

            try:
                get_mlst_files(url, str(scheme_dir), client_key, client_secret,
                               session_token, session_secret, scheme, verbose=verbose,
                               no_auth=no_auth, session=session,
                               show_progress=not parallel)
                return ('ok', scheme, database)
            except Exception as e:
                return ('error', scheme, database, e)

        if parallel:
            info(f"Downloading {len(scheme_entries)} schemes using {threads} threads...")

        with ThreadPoolExecutor(max_workers=threads) as executor:
            future_to_entry = {
                executor.submit(_download_scheme, entry): entry
                for entry in scheme_entries
            }
            for future in tqdm(as_completed(future_to_entry),
                              total=len(scheme_entries),
                              desc="Downloading MLST schemes", unit="scheme"):
                result = future.result()
                status = result[0]
                scheme = result[1]
                database = result[2]

                if status == 'ok':
                    success(f"Successfully downloaded scheme: {scheme}")
                    download_success = True
                elif status == 'skip_cred':
                    error(f"No credentials available for {database}, skipping {scheme}")
                    skipped_schemes.append((scheme, database))
                elif status == 'error':
                    exc = result[3]
                    if '401 Client Error' in str(exc) or '403 Client Error' in str(exc):
                        if no_auth:
                            error(f"Scheme [{scheme}] requires authentication and cannot be downloaded with --no-auth.")
                            info(f"Run without --no-auth to download authenticated schemes.")
                        else:
                            error(f"Authentication error for [{scheme}]: {exc}")
                            error(f"Please make sure you have registered to [{scheme}] scheme within the {database} database.")
                    else:
                        error(f"Error downloading scheme [{scheme}]: {exc}")
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