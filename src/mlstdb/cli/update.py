import click
from pathlib import Path
from mlstdb.core.auth import get_client_credentials, retrieve_session_token
from mlstdb.core.download import get_mlst_files, create_blast_db
from mlstdb.core.config import check_dir
from mlstdb.utils import error, success, info
from tqdm import tqdm
import sys

@click.command()
@click.help_option('-h', '--help')
@click.option('--input', '-i', required=True, 
              help='Path to mlst_schemes_<db>.tab containing MLST scheme URLs')
@click.option('--directory', '-d', default='pubmlst',
              help='Directory to save the downloaded MLST schemes (default: pubmlst)')
@click.option('--blast-directory', '-b',
              help='Directory for BLAST database (default: blast)')
@click.option('--verbose', '-v', is_flag=True,
              help='Enable verbose logging for debugging')

def update(input: str, directory: str, blast_directory: str, verbose: bool):
    """
    Update MLST schemes and create BLAST database.

    Downloads MLST schemes from the specified input file and creates a BLAST database
    from the downloaded sequences. Authentication tokens should be set up using fetch.py.
    """
    try:
        # Read the input file
        with open(input, 'r') as f:
            # Skip header
            header = next(f)
            lines = f.readlines()

        check_dir(directory)
        
        # Track credential issues
        auth_errors = False
        download_success = False

        # Process each scheme
        for line in tqdm(lines, desc="Downloading MLST schemes", unit="scheme"):
            parts = line.strip().split('\t')
            if len(parts) != 5:
                error(f"Skipping invalid line: {line}")
                continue

            database, species, scheme_desc, scheme, url = parts
            
            try:
                # Get credentials for the specific database
                try:
                    client_key, client_secret = get_client_credentials(database.lower())
                    session_token, session_secret = retrieve_session_token(database.lower())
                except ValueError as ve:
                    error(f"Error with credentials for {database}: {ve}")
                    auth_errors = True
                    continue

                if not session_token or not session_secret:
                    error(f"No valid session token found for {database}.")
                    auth_errors = True
                    continue

                scheme_dir = Path(directory) / scheme
                check_dir(str(scheme_dir))

                try:
                    get_mlst_files(url, str(scheme_dir), client_key, client_secret,
                                session_token, session_secret, scheme,
                                verbose=verbose)
                    success(f"Successfully downloaded scheme: {scheme}")
                    download_success = True
                except Exception as e:
                    # Check for 401/403 errors specifically
                    if '401 Client Error: Unauthorized' in str(e) or '403 Client Error: Forbidden' in str(e):
                        error(f"Authentication error for {scheme}: {e}")
                        auth_errors = True
                    else:
                        error(f"Error downloading scheme {scheme}: {e}")
                    continue
            
            except Exception as e:
                error(f"Error downloading scheme {scheme}: {e}")
                continue
        
        # Exit early if we have authentication errors
        if auth_errors:
            error("\nAuthentication errors occurred during downloads.")
            info("\nTo fix authentication issues:")
            info("1. Run 'mlstdb fetch' to refresh or setup your credentials")
            info("2. Then run this command again")
            sys.exit(1)

        # Check if we have any schemes downloaded
        scheme_dirs = [d for d in Path(directory).iterdir() if d.is_dir()]
        if not scheme_dirs:
            error("\nNo schemes were successfully downloaded. BLAST database creation skipped.")
            sys.exit(1)

        # Create BLAST database after all schemes are downloaded
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