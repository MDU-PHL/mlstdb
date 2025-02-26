#!/usr/bin/env python3
import os
import requests
import re
import click
from pathlib import Path
from tqdm import tqdm
from rauth import OAuth1Session, OAuth1Service
import sys
import logging
import configparser
from typing import List, Dict, Set, Tuple, Optional

# Constants
BASE_WEB = {
    "pubmlst": "https://pubmlst.org/bigsdb",
    "pasteur": "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl",
}

BASE_API = {
    "pubmlst": "https://rest.pubmlst.org",
    "pasteur": "https://bigsdb.pasteur.fr/api",
}

DB_MAPPING = {
    "pubmlst": "pubmlst_neisseria_seqdef",
    "pasteur": "pubmlst_diphtheria_isolates"
}

# Color formatting using click
def error(msg: str) -> None:
    click.secho(f"Error: {msg}", fg="red", err=True)

def success(msg: str) -> None:
    click.secho(msg, fg="green")

def info(msg: str) -> None:
    click.secho(msg, fg="blue")

def get_config_dir() -> Path:
    """Create and return the configuration directory."""
    config_dir = Path.home() / ".config" / "mlstdb"
    if not config_dir.exists():
        config_dir.mkdir(parents=True, mode=0o700)
    return config_dir

def setup_client_credentials(site: str) -> Tuple[str, str]:
    """Setup and save client credentials."""
    config = configparser.ConfigParser(interpolation=None)
    file_path = get_config_dir() / "client_credentials"
    
    if file_path.exists():
        config.read(file_path)
    
    info("\nPlease enter your client credentials:")
    client_id = click.prompt("Client ID", type=str).strip()
    while len(client_id) != 24:
        error("Client IDs must be exactly 24 characters long")
        client_id = click.prompt("Client ID", type=str).strip()
    
    client_secret = click.prompt("Client Secret", type=str).strip()
    while len(client_secret) != 42:
        error("Client secrets must be exactly 42 characters long")
        client_secret = click.prompt("Client Secret", type=str).strip()

    config[site] = {"client_id": client_id, "client_secret": client_secret}
    
    with open(file_path, "w") as configfile:
        config.write(configfile)
    success(f"\nClient credentials saved to {file_path}")
    return client_id, client_secret

def register_tokens(db: str):
    """Setup authentication tokens by registering with the service."""
    info(f"\nNo tokens found for {db}. Starting registration process...")
    
    # Setup client credentials
    client_id, client_secret = setup_client_credentials(db)
    
    # Initialize OAuth service
    service = OAuth1Service(
        name="BIGSdb_downloader",
        consumer_key=client_id,
        consumer_secret=client_secret,
        request_token_url=f"{BASE_API[db]}/db/{DB_MAPPING[db]}/oauth/get_request_token",
        access_token_url=f"{BASE_API[db]}/db/{DB_MAPPING[db]}/oauth/get_access_token",
        base_url=BASE_API[db],
    )
    
    # Get request token
    info("\nRequesting temporary token...")
    r = service.get_raw_request_token(
        params={"oauth_callback": "oob"},
        headers={"User-Agent": "BIGSdb downloader"}
    )
    if r.status_code != 200:
        error(f"Failed to get request token: {r.json()['message']}")
        sys.exit(1)
    
    request_token = r.json()["oauth_token"]
    request_secret = r.json()["oauth_token_secret"]
    success("Temporary token received")
    
    # Get access token
    click.secho("\nAuthorization Required", fg="yellow", bold=True)
    info(
        "\nPlease open this URL in your browser:\n"
        f"{BASE_WEB[db]}?db={DB_MAPPING[db]}&page=authorizeClient&oauth_token={request_token}"
    )
    
    verifier = click.prompt("\nEnter the verification code from the website", type=str)
    
    info("\nRequesting access token...")
    r = service.get_raw_access_token(
        request_token,
        request_secret,
        params={"oauth_verifier": verifier},
        headers={"User-Agent": "BIGSdb downloader"},
    )
    
    if r.status_code != 200:
        error(f"Failed to get access token: {r.json()['message']}")
        sys.exit(1)
        
    access_token = r.json()["oauth_token"]
    access_secret = r.json()["oauth_token_secret"]
    
    # Save access token
    config = configparser.ConfigParser(interpolation=None)
    file_path = get_config_dir() / "access_tokens"
    if file_path.exists():
        config.read(file_path)
    config[db] = {"token": access_token, "secret": access_secret}
    with open(file_path, "w") as configfile:
        config.write(configfile)
    success(f"\nAccess token saved to {file_path}")

    # Get session token
    info("\nRequesting session token...")
    url = f"{BASE_API[db]}/db/{DB_MAPPING[db]}/oauth/get_session_token"
    
    session = OAuth1Session(
        service.consumer_key,
        service.consumer_secret,
        access_token=access_token,
        access_token_secret=access_secret
    )
    
    r = session.get(url, headers={"User-Agent": "BIGSdb downloader"})
    
    if r.status_code != 200:
        error(f"Failed to get session token: {r.json()['message']}")
        sys.exit(1)
        
    token = r.json()["oauth_token"]
    secret = r.json()["oauth_token_secret"]
    
    # Save session token
    config = configparser.ConfigParser(interpolation=None)
    file_path = get_config_dir() / "session_tokens"
    if file_path.exists():
        config.read(file_path)
    config[db] = {"token": token, "secret": secret}
    with open(file_path, "w") as configfile:
        config.write(configfile)
    
    success(f"\nSession token saved to {file_path}")
    return token, secret

@click.command()
@click.option('--db', '-d', type=click.Choice(['pubmlst', 'pasteur']), 
              help='Database to use (pubmlst or pasteur)')
@click.option('--exclude', '-e', default='cgMLST', 
              help='Scheme name must not include provided term (default: cgMLST)')
@click.option('--match', '-m', default='MLST', 
              help='Scheme name must include provided term (default: MLST)')
@click.option('--scheme-uris', '-s', default='scheme_uris.tab',
              help='Path to scheme_uris.tab file for scheme sanitisation')
@click.option('--filter', '-f',
              help='Filter species or schemes using a wildcard pattern during scheme sanitisation')
@click.option('--resume', '-r', is_flag=True, 
              help='Resume processing from where it stopped')
@click.option('--verbose', '-v', is_flag=True, 
              help='Enable verbose logging for debugging')
def main(db, exclude, match, scheme_uris, filter, resume, verbose):  # Added scheme_uris parameter
    """BIGSdb Scheme Fetcher Tool
    
    This tool downloads MLST scheme information from BIGSdb databases.
    It will automatically handle authentication and save the results.
    """
    # If db is not provided, prompt for it
    if not db:
        db = click.prompt(
            "Which database would you like to use?",
            type=click.Choice(['pubmlst', 'pasteur']),
            default='pubmlst'
        )

    config_dir = get_config_dir()
    client_creds_file = config_dir / "client_credentials"
    session_tokens_file = config_dir / "session_tokens"

    # Check if credentials exist, if not setup
    if not client_creds_file.exists() or not session_tokens_file.exists():
        register_tokens(db)

    # Get credentials
    import configparser
    config = configparser.ConfigParser(interpolation=None)
    
    # Read client credentials
    config.read(client_creds_file)
    if not config.has_section(db):
        error(f"No client credentials found for {db}")
        register_tokens(db)
        config.read(client_creds_file)
    
    client_key = config[db]["client_id"]
    client_secret = config[db]["client_secret"]

    # Read session tokens
    config.read(session_tokens_file)
    if not config.has_section(db):
        error(f"No session token found for {db}")
        register_tokens(db)
        config.read(session_tokens_file)
    
    session_token = config[db]["token"]
    session_secret = config[db]["secret"]

    output_file = f"mlst_schemes_{db}.txt"
    processed_file = f"processed_dbs_{db}.txt"

    try:
        if not resume:
            clear_file(output_file)
            clear_file(processed_file)
        
        processed_dbs = load_processed_databases(processed_file)
        
        base_uri = BASE_API[db]
        resources = fetch_resources(base_uri, client_key, client_secret, 
                                  session_token, session_secret, verbose)
        
        success_count = 0
        total_dbs = sum(1 for r in resources if 'databases' in r 
                       for _ in r['databases'])
        
        for resource in tqdm(resources, desc="Processing resources"):
            if 'databases' in resource:
                for database in resource['databases']:
                    if database['description'] in processed_dbs:
                        if verbose:
                            info(f"Skipping already processed database: {database['description']}")
                        success_count += 1
                        continue
                    
                    try:
                        get_matching_schemes(database, match, exclude, 
                                           client_key, client_secret,
                                           session_token, session_secret,
                                           output_file, processed_file, verbose)
                        success_count += 1
                    except Exception as e:
                        error(f"Error processing {database['description']}: {e}")
                        continue

        # Delete processed_file if all databases were processed successfully
        if success_count == total_dbs:
            try:
                os.remove(processed_file)
                if verbose:
                    info(f"Removed progress tracking file: {processed_file}")
            except OSError as e:
                error(f"Error removing progress file: {e}")

        # After successful fetch, perform scheme sanitisation
        if Path(scheme_uris).exists():
            sanitise_output(output_file, scheme_uris, filter, verbose)
        else:
            error(f"Scheme URIs file not found: {scheme_uris}")
            error("Skipping scheme sanitisation step")    
                    
    except Exception as e:
        error(f"An error occurred: {e}")
        info("Progress file kept for resume capability")
        sys.exit(1)

def clear_file(file_path):
    """Clear the contents of a file or create it if it doesn't exist."""
    with open(file_path, 'w') as f:
        # Add headers only for mlst_schemes files
        if file_path.startswith('mlst_schemes_'):
            f.write("database\tspecies\tscheme_description\tURI\n")

        
def fetch_resources(base_uri, client_key, client_secret, session_token, session_secret, verbose=False):
    if verbose:
        print(f"Fetching resources from {base_uri}")
    return fetch_json(base_uri, client_key, client_secret, session_token, session_secret, verbose)

def load_processed_databases(file_path):
    if not os.path.exists(file_path):
        return set()
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f)

def save_processed_database(file_path, db_description):
    with open(file_path, 'a') as f:
        f.write(f"{db_description}\n")

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
        
        # Get database type from output filename
        db_type = output_file.split('_')[-1].replace('.txt', '')
        
        if schemes and 'schemes' in schemes:
            matching_schemes = []
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


def remove_db_credentials(config_dir: Path, db: str) -> None:
    """Remove credentials for specific database while preserving others."""
    for file_name in ["client_credentials", "session_tokens", "access_tokens"]:
        file_path = config_dir / file_name
        if file_path.exists():
            config = configparser.ConfigParser(interpolation=None)
            config.read(file_path)
            if db in config:
                config.remove_section(db)
                with open(file_path, 'w') as f:
                    config.write(f)
                success(f"Removed {db} credentials from {file_name}")


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
        
        # Handle 401 Unauthorized error - try once to refresh token
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
                    info("Please run the script again to generate new credentials.")
                    sys.exit(1)
                except Exception as del_error:
                    error(f"Failed to delete credentials: {del_error}")
                    sys.exit(1)
            
            error("Exiting. Please fix credentials and try again")
            sys.exit(1)

        raise

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

    info(f"Sanitising output file: {output_file}")
    
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
      
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()