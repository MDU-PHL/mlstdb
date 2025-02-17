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
              help='Path to mlst_schemes_<db>.txt containing MLST scheme URLs')
@click.option('--directory', '-d', default='pubmlst',
              help='Directory to save the downloaded MLST schemes (default: pubmlst)')
@click.option('--blast-directory', '-b',
              help='Directory for BLAST database (default: relative to directory as ../blast)')
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

        # Process each scheme
        for line in tqdm(lines, desc="Downloading MLST schemes", unit="scheme"):
            parts = line.strip().split('\t')
            if len(parts) != 5:
                error(f"Skipping invalid line: {line}")
                continue

            database, species, scheme_desc, scheme, url = parts
            
            try:
                # Get credentials for the specific database
                client_key, client_secret = get_client_credentials(database.lower())
                session_token, session_secret = retrieve_session_token(database.lower())

                if not session_token or not session_secret:
                    error(f"No valid session token found for {database}. Please run fetch.py first to set up authentication.")
                    continue

                # scheme_dir = Path(directory) / sanitise_name(scheme)
                scheme_dir = Path(directory) / scheme
                check_dir(str(scheme_dir))

                get_mlst_files(url, str(scheme_dir), client_key, client_secret,
                             session_token, session_secret, scheme,
                             verbose=verbose)
                success(f"Successfully downloaded scheme: {scheme}")

            except Exception as e:
                error(f"Error downloading scheme {scheme}: {e}")
                if verbose:
                    import traceback
                    error(traceback.format_exc())
                continue

        # Create BLAST database after all schemes are downloaded
        info("\nCreating BLAST database from downloaded MLST schemes...")
        create_blast_db(directory, blast_directory, verbose)
        success("Update completed successfully!")

    except Exception as e:
        error(f"An error occurred: {e}")
        if verbose:
            import traceback
            error(traceback.format_exc())
        sys.exit(1)