import os
import requests
import click
import configparser
from pathlib import Path
from tqdm import tqdm
from rauth import OAuth1Session
import logging
import sys
from typing import Tuple, List
import subprocess
import glob


BASE_API = {
    "pubmlst": "https://rest.pubmlst.org",
    "pasteur": "https://bigsdb.pasteur.fr/api",
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

def get_client_credentials(key_name: str) -> Tuple[str, str]:
    """Get OAuth client credentials from config file."""
    config = configparser.ConfigParser(interpolation=None)
    file_path = get_config_dir() / "client_credentials"
    
    if file_path.is_file():
        config.read(file_path)
        if config.has_section(key_name):
            return (config[key_name]["client_id"], 
                   config[key_name]["client_secret"])
    
    raise ValueError(f"Client credentials not found for {key_name}")

def retrieve_session_token(key_name: str) -> Tuple[str, str]:
    """Get OAuth session token from config file."""
    config = configparser.ConfigParser(interpolation=None)
    file_path = get_config_dir() / "session_tokens"
    
    if file_path.is_file():
        config.read(file_path)
        if config.has_section(key_name):
            return (config[key_name]["token"], 
                   config[key_name]["secret"])
    
    return None, None

def check_dir(directory: str) -> None:
    """Ensure the directory exists and is writable."""
    path = Path(directory)
    if not path.exists():
        path.mkdir(parents=True)
    if not (path.is_dir() and os.access(directory, os.W_OK)):
        raise PermissionError(f"Cannot write to directory: {directory}")

def sanitise_name(name: str) -> str:
    """Sanitise directory or file names by replacing invalid characters."""
    return name.replace('/', '_').replace('\\', '_').replace(':', '_')

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
        profiles_file_path = os.path.join(directory, f"{sanitise_name(scheme_name)}.txt")
        with open(profiles_file_path, 'w') as f:
            f.write(profiles.text)
            
    except requests.exceptions.HTTPError as e:
        if e.response.status_code in [401, 403]:
            error("\nAuthentication failed!")
            info("\nTo fix authentication issues:")
            info("Run 'mlstdb fetch' to refresh your credentials for the database you are trying to access.")
            info("Then try running this script again.")
            # info("2. Try running this script again")
            sys.exit(1)
        elif e.response.status_code == 404:
            error(f"Resource not found at URL: {url}")
        raise

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

@click.command()
@click.option('--input', '-i', required=True, 
              help='Path to mlst_schemes_<db>.txt containing MLST scheme URLs')
@click.option('--directory', '-d', default='pubmlst',
              help='Directory to save the downloaded MLST schemes (default: pubmlst)')
@click.option('--blast-directory', '-b',
              help='Directory for BLAST database (default: relative to directory as ../blast)')
@click.option('--verbose', '-v', is_flag=True,
              help='Enable verbose logging for debugging')
def main(input: str, directory: str, blast_directory: str, verbose: bool):
    """
    BIGSdb MLST Scheme Downloader Tool

    This tool downloads MLST schemes and their associated files from BIGSdb databases,
    then creates a BLAST database from the downloaded sequences.
    It requires authentication tokens which can be set up using the fetch.py script.
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

                scheme_dir = Path(directory) / sanitise_name(scheme)
                check_dir(scheme_dir)

                get_mlst_files(url, scheme_dir, client_key, client_secret,
                             session_token, session_secret, scheme,
                             verbose=verbose)
                success(f"Successfully downloaded scheme: {scheme}")

            except Exception as e:
                error(f"Error downloading scheme {scheme}: {e}")
                continue

        # Create BLAST database after all schemes are downloaded
        info("\nCreating BLAST database from downloaded MLST schemes...")
        create_blast_db(directory, blast_directory, verbose)

    except Exception as e:
        error(f"An error occurred: {e}")
        sys.exit(1)
        
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()