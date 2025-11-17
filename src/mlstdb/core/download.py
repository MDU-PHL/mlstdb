import os
import re
import sys
import click
import requests
import subprocess
import configparser
from pathlib import Path
from tqdm import tqdm
from rauth import OAuth1Session, OAuth1Service
from mlstdb.core.config import get_config_dir, BASE_API, DB_MAPPING
from mlstdb.core.auth import remove_db_credentials, register_tokens
from mlstdb.utils import error, success, info

def fetch_json(url, client_key, client_secret, session_token, session_secret, verbose=False):
    """Fetch JSON from URL with OAuth authentication and session token refresh."""
    if verbose:
        print(f"Fetching JSON from {url}")
    
    # Initialize session with current token
    session = OAuth1Session(
        consumer_key=client_key,
        consumer_secret=client_secret,
        access_token=session_token,
        access_token_secret=session_secret,
    )
    session.headers.update({"User-Agent": "BIGSdb downloader"})

    try:
        response = session.get(url)
        if verbose:
            print(f"Response code: {response.status_code}, URL: {url}")
        
        if response.status_code == 404:
            print(f"Resource not found at URL: {url}")
            return None
        
        # Handle 401 Unauthorised error - try once to refresh token
        if response.status_code == 401:
            info("Invalid session token. Requesting new one...")
            
            # Determine which database we're working with
            if url.startswith(BASE_API['pubmlst']):
                db = 'pubmlst'
            elif url.startswith(BASE_API['pasteur']):
                db = 'pasteur'
            else:
                raise ValueError(f"Unable to determine database from URL: {url}")
            
            # Get new session token using existing credentials
            config = configparser.ConfigParser(interpolation=None)
            access_tokens_file = get_config_dir() / "access_tokens"
            
            # Read access tokens
            config.read(access_tokens_file)
            access_token = config[db]["token"]
            access_secret = config[db]["secret"]
            
            # Initialize OAuth service
            service = OAuth1Service(
                name="BIGSdb_downloader",
                consumer_key=client_key,
                consumer_secret=client_secret,
                request_token_url=f"{BASE_API[db]}/db/{DB_MAPPING[db]}/oauth/get_request_token",
                access_token_url=f"{BASE_API[db]}/db/{DB_MAPPING[db]}/oauth/get_access_token",
                base_url=BASE_API[db],
            )
            
            # Get new session token
            url_session = f"{BASE_API[db]}/db/{DB_MAPPING[db]}/oauth/get_session_token"
            session_request = OAuth1Session(
                client_key,
                client_secret,
                access_token=access_token,
                access_token_secret=access_secret,
            )
            session_request.headers.update({"User-Agent": "BIGSdb downloader"})
            
            r = session_request.get(url_session)
            if r.status_code == 200:
                new_token = r.json()["oauth_token"]
                new_secret = r.json()["oauth_token_secret"]
                
                # Save new session token
                config = configparser.ConfigParser(interpolation=None)
                session_tokens_file = get_config_dir() / "session_tokens"
                if session_tokens_file.exists():
                    config.read(session_tokens_file)
                config[db] = {"token": new_token, "secret": new_secret}
                with open(session_tokens_file, "w") as configfile:
                    config.write(configfile)
                
                if verbose:
                    success("New session token obtained and saved")
                
                info("\nSession token has been refreshed. Please run the command again.")
                sys.exit(0)  # Exit cleanly after token refresh
            else:
                # If we can't get a new session token, raise the original 401 error
                response.raise_for_status()
        
        response.raise_for_status()
        return response.json()
        
    except requests.exceptions.HTTPError as e:
        if e.response.status_code in [401, 403]:  # Only handle other 401/403 cases
            config_dir = get_config_dir()
            
            error(f"\nAuthentication Failed! (Status code: {e.response.status_code})")
            if e.response.status_code == 401:
                info("This usually means your tokens have expired. Please check your credentials.")
            else:
                info("This usually means you lack permissions for this resource. Please check your credentials.")
                
            # Show credential locations
            info("\nYour credentials are stored in:")
            for cred_file in ["client_credentials", "session_tokens", "access_tokens"]:
                info(f"- {config_dir}/{cred_file}")

            # Inform user about next steps
            info("\nTo fix authentication issues:")
            info("1. You need to manually delete your credentials from the files above")
            info("2. Run the script again to generate new credentials")
            
            # Offer to delete credentials
            if click.confirm("\nWould you like to delete credentials for this database?", default=False):
                try:
                    # Use the db parameter from the main function
                    if url.startswith(BASE_API['pubmlst']):
                        db = 'pubmlst'
                    elif url.startswith(BASE_API['pasteur']):
                        db = 'pasteur'
                    remove_db_credentials(config_dir, db)
                    info("\nCredentials deleted successfully.")
                    info("Please run the command again to generate new credentials.")
                    sys.exit(1)
                except Exception as del_error:
                    error(f"Failed to delete credentials: {del_error}")
                    sys.exit(1)
            
            error("Exiting. Please fix credentials and try again")
            sys.exit(1)

        raise


def get_mlst_files(url: str, directory: str, client_key: str, client_secret: str, 
                   session_token: str, session_secret: str, scheme_name: str, 
                   verbose: bool = False) -> None:
    """Download MLST data and save them in the given directory."""
    session = OAuth1Session(
        consumer_key=client_key,
        consumer_secret=client_secret,
        access_token=session_token,
        access_token_secret=session_secret,
    )
    
    session.headers.update({"User-Agent": "BIGSdb downloader"})

    if verbose:
        info(f"Fetching MLST scheme from {url}...")

    try:
        response = session.get(url)    
        response.raise_for_status()
        mlst_scheme = response.json()
        
        if verbose:
            info(f"Retrieved MLST scheme: {mlst_scheme}")

        # Modified version handling for null values
        db_version = mlst_scheme.get('last_added', 'Not found')
        if db_version is None:
            db_version = 'No version information available'
        info(f"Database version: {db_version}")

        # Save database version to a file
        db_version_path = os.path.join(directory, 'database_version.txt')
        with open(db_version_path, 'w') as version_file:
            version_file.write(db_version)

        # Download loci with progress bar
        for loci in tqdm(mlst_scheme['loci'], desc="Downloading loci", unit="locus"):
            name = loci.split('/')[-1]
            loci_fasta = session.get(loci + '/alleles_fasta')
            loci_fasta.raise_for_status()
            loci_file_name = os.path.join(directory, name + '.tfa')
            with open(loci_file_name, 'wb') as f:
                f.write(loci_fasta.content)

        # Download profiles CSV
        profiles_url = url + '/profiles_csv'
        profiles = session.get(profiles_url)
        profiles.raise_for_status()
        profiles_file_path = os.path.join(directory, f"{scheme_name}.txt")
        with open(profiles_file_path, 'w') as f:
            f.write(profiles.text)
            
        # Handle expired token (401)
        if response.status_code == 401:
            if verbose:
                info("Session token expired, attempting to refresh...")
            # Get database type from URL    
            db = 'pubmlst' if 'pubmlst.org' in url else 'pasteur'
            
            # Register new tokens
            new_token, new_secret = register_tokens(db)
            
            # Retry with new token
            session = OAuth1Session(
                consumer_key=client_key,
                consumer_secret=client_secret,
                access_token=new_token,
                access_token_secret=new_secret,
            )
            session.headers.update({"User-Agent": "BIGSdb downloader"})
            response = session.get(url)
                
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 403:
            error("\nAuthentication failed - permission denied!")
            info("\nTo fix permission issues:")
            info("Run 'mlstdb fetch' to refresh your credentials for the database you are trying to access.")
            info("Then try running this script again.")
            sys.exit(1)
        elif e.response.status_code == 404:
            error(f"Resource not found at URL: {url}")
        raise

def fetch_resources(base_uri, client_key, client_secret, session_token, session_secret, verbose=False):
    if verbose:
        print(f"Fetching resources from {base_uri}")
    return fetch_json(base_uri, client_key, client_secret, session_token, session_secret, verbose)

def clear_file(file_path):
    """Clear the contents of a file or create it if it doesn't exist."""
    with open(file_path, 'w') as f:
        # Add headers only for mlst_schemes files
        if file_path.startswith('mlst_schemes_'):
            f.write("database\tspecies\tscheme_description\tURI\n")

def load_processed_databases(file_path):
    if not os.path.exists(file_path):
        return set()
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f)

def save_processed_database(file_path, db_description):
    with open(file_path, 'a') as f:
        f.write(f"{db_description}\n")


def sanitise_species(species_column: str) -> str:
    """Clean up species names."""
    return species_column.replace("sequence/profile definitions", "").strip()

def extract_scheme_from_url(url: str) -> str:
    """Extract scheme name from URL pattern."""
    # Handle Pasteur URLs
    if 'pasteur.fr' in url:
        match = re.search(r'pubmlst_([^_]+)_seqdef/schemes/(\d+)', url)
        if match:
            scheme = match.group(1)
            scheme_id = match.group(2)
            return f"{scheme}_{scheme_id}" if scheme_id != "1" else scheme
    
    # Handle PubMLST URLs
    else:
        match = re.search(r'pubmlst_([^_]+)_seqdef/schemes/(\d+)', url)
        if match:
            scheme = match.group(1)
            scheme_id = match.group(2)
            return f"{scheme}_{scheme_id}" if scheme_id != "1" else scheme
    
    return "missing"

def load_scheme_uris(file_path: str) -> dict:
    """Load scheme URIs from tab-delimited file."""
    schemes = {}
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            scheme, uri = line.strip().split('\t')
            schemes[uri] = scheme
    return schemes

def sanitise_output(output_file: str, scheme_uris_file: str, filter_pattern: str, verbose: bool) -> None:
    """Sanitise the output file using the scheme URIs mapping."""
    if not Path(output_file).exists():
        error(f"Output file not found: {output_file}")
        return

    scheme_uris = load_scheme_uris(scheme_uris_file)
    sanitised_data = []
    existing_schemes = set()

    info(f"Sanitising schemes of the output file: {output_file}")
    
    # Read and process the file
    with open(output_file, 'r') as infile:
        header = next(infile)  # Skip header
        for line in infile:
            columns = line.strip().split('\t')
            if len(columns) != 4:  # Now expecting 4 columns
                error(f"Skipping malformed line: {line}")
                continue

            database = columns[0]  # database type (pubmlst/pasteur)
            species = sanitise_species(columns[1])
            scheme_desc = columns[2]
            uri = columns[3]

            if filter_pattern and not (
                re.search(filter_pattern, species, re.IGNORECASE) or
                re.search(filter_pattern, uri, re.IGNORECASE)
            ):
                if verbose:
                    info(f"Skipping entry due to filter mismatch: {line.strip()}")
                continue

            # First try to get scheme from mapping file
            scheme = scheme_uris.get(uri)
            existing_schemes.add(scheme if scheme else "missing")
            sanitised_data.append((database, species, scheme_desc, scheme, uri))

    # Handle missing schemes
    if any(entry[3] is None for entry in sanitised_data):
        error("\nThe following URIs have missing schemes:")
        for entry in sanitised_data:
            if entry[3] is None:
                click.echo(entry[4])
        
        user_choice = click.prompt(
            "\nDo you want to set missing schemes as 'missing' or auto-generate them?",
            type=click.Choice(['missing', 'auto'], case_sensitive=False)
        )

        if user_choice == "auto":
            for idx, (database, species, scheme_desc, scheme, uri) in enumerate(sanitised_data):
                if scheme is None:
                    auto_scheme = extract_scheme_from_url(uri)
                    if verbose:
                        info(f"Auto-generated scheme: {auto_scheme} for URI: {uri}")
                    sanitised_data[idx] = (database, species, scheme_desc, auto_scheme, uri)
        else:  # missing
            sanitised_data = [
                (database, species, scheme_desc, scheme if scheme else "missing", uri)
                for database, species, scheme_desc, scheme, uri in sanitised_data
            ]

    # Write sanitised output back to the same file
    with open(output_file, 'w') as outfile:
        outfile.write("database\tspecies\tscheme_description\tscheme\tURI\n")
        for entry in sanitised_data:
            outfile.write('\t'.join(str(x) for x in entry) + '\n')
    
    success(f"Scheme sanitisation complete! Results updated in {output_file}")


def get_matching_schemes(db, match, exclude, client_key, client_secret, 
                        session_token, session_secret, output_file, processed_file, verbose=False):
    """Get matching schemes from the database."""
    # Check if this is a sequence definition database
    is_seqdef = ('seqdef' in db['name'].lower() or 
                 'definitions' in db.get('description', '').lower())
    
    if is_seqdef:
        try:
            db_attributes = fetch_json(db['href'], client_key, client_secret, 
                                     session_token, session_secret, verbose=verbose)
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 401:
                print(f"The token does not allow access to the `{db['description']}` database. "
                      f"So, the `{db['description']}` database will be skipped.\n"
                      "To download the data, please ensure your account has access to this database.")
                save_processed_database(processed_file, db['description'])
                return
            elif e.response.status_code == 404:
                print(f"The resource `{db['href']}` for the `{db['description']}` database was not found (404). "
                      "Skipping this database.")
                save_processed_database(processed_file, db['description'])
                return
            else:
                raise
        
        if not db_attributes or 'schemes' not in db_attributes:
            save_processed_database(processed_file, db['description'])
            return
        
        schemes = fetch_json(db_attributes['schemes'], client_key, client_secret, 
                           session_token, session_secret, verbose=verbose)
             
        if schemes and 'schemes' in schemes:
            matching_schemes = []
            
            # Determine database type from URL instead of filename
            db_type = "pasteur" if "pasteur.fr" in db['href'] else "pubmlst"

            for scheme in schemes['schemes']:
                if match and not re.search(match, scheme['description'], flags=0):
                    continue
                if exclude and re.search(exclude, scheme['description'], flags=0):
                    continue
                # Add database type as first column
                matching_schemes.append(f"{db_type}\t{db['description']}\t{scheme['description']}\t{scheme['scheme']}\n")
            
            if matching_schemes:  # Only write to file if we found matching schemes
                with open(output_file, 'a') as f:
                    f.writelines(matching_schemes)
        
        save_processed_database(processed_file, db['description'])

def create_blast_db(input_dir: str, blast_directory: str, verbose: bool = False) -> None:
    """Create BLAST database from MLST schemes."""
    input_path = Path(input_dir)
    if not input_path.is_dir():
        error(f"Input directory {input_dir} does not exist.")
        return

    # Set default output directory relative to input
    blast_path = Path(blast_directory) if blast_directory else input_path.parent / "blast"
    blast_path.mkdir(parents=True, exist_ok=True)

    blast_file = blast_path / "mlst.fa"
    if blast_file.exists():
        blast_file.unlink()

    info(f"Creating BLAST database from {input_dir}")
    info(f"Output directory: {blast_path}")

    # Get all scheme directories
    scheme_dirs = [d for d in input_path.iterdir() if d.is_dir()]
    total_schemes = len(scheme_dirs)

    if total_schemes == 0:
        error("No scheme directories found.")
        return

    with open(blast_file, 'w') as outfile:
        for scheme_dir in tqdm(scheme_dirs, desc="Processing schemes", unit="scheme"):
            scheme_name = scheme_dir.name
            if verbose:
                info(f"\nProcessing scheme: {scheme_name}")

            # Find all sequence files
            sequence_files = []
            for ext in ['.tfa', '.fasta', '.fa', '.fas']:
                sequence_files.extend(scheme_dir.glob(f"*{ext}"))

            if not sequence_files:
                if verbose:
                    info(f"No sequence files found in {scheme_dir}")
                continue

            # Process each sequence file
            for seq_file in sequence_files:
                with open(seq_file) as f:
                    for line in f:
                        if line.startswith('>'):
                            # Modify header to include scheme name
                            outfile.write(f">{scheme_name}.{line[1:]}")
                        else:
                            outfile.write(line)

    if blast_file.stat().st_size == 0:
        error("No sequences found to create BLAST database.")
        return

    # Create BLAST database
    info("\nCreating BLAST database...")
    try:
        subprocess.run([
            'makeblastdb',
            '-hash_index',
            '-in', str(blast_file),
            '-dbtype', 'nucl',
            '-title', 'MLST',
            '-parse_seqids'
        ], check=True, capture_output=True, text=True)
        success(f"BLAST database created successfully: {blast_file}")
    except subprocess.CalledProcessError as e:
        error(f"Failed to create BLAST database: {e.stderr}")
    except FileNotFoundError:
        error("makeblastdb command not found. Please ensure BLAST+ is installed.")