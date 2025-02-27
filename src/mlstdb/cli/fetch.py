import click
from pathlib import Path
import importlib.resources
from tqdm import tqdm
import configparser
import sys
import os
from mlstdb.core.auth import register_tokens, setup_client_credentials, remove_db_credentials
from mlstdb.core.download import (fetch_resources, get_matching_schemes, 
                                sanitise_output, clear_file, 
                                load_processed_databases,save_processed_database, load_scheme_uris)
from mlstdb.core.config import get_config_dir, BASE_API
from mlstdb.utils import error, success, info


@click.command()
@click.help_option('-h', '--help')
@click.option('--db', '-d', type=click.Choice(['pubmlst', 'pasteur']), 
              help='Database to use (pubmlst or pasteur)')
@click.option('--exclude', '-e', default='cgMLST', 
              help='Scheme name must not include provided term (default: cgMLST)')
@click.option('--match', '-m', default='MLST', 
              help='Scheme name must include provided term (default: MLST)')
@click.option('--scheme-uris', '-s',
              help='Optional: Path to custom scheme_uris.tab file')
@click.option('--filter', '-f',
              help='Filter species or schemes using a wildcard pattern')
@click.option('--resume', '-r', is_flag=True, 
              help='Resume processing from where it stopped')
@click.option('--verbose', '-v', is_flag=True, 
              help='Enable verbose logging for debugging')
def fetch(db, exclude, match, scheme_uris, filter, resume, verbose):
    """BIGSdb Scheme Fetcher Tool
    
    This tool downloads MLST scheme information from BIGSdb databases.
    It will automatically handle authentication and save the results.
    """
    
    try:
        # If scheme_uris is not provided, use the package data
        if scheme_uris is None:
            with importlib.resources.path('mlstdb.data', 'scheme_uris.tab') as default_path:
                scheme_uris = str(default_path)
        
        # If db is not provided, prompt for it
        if not db:
            db = click.prompt(
                "Which database would you like to use?",
                type=click.Choice(['pubmlst', 'pasteur']),
                default='pubmlst'
            )
        
        # Get client credentials
        config_dir = get_config_dir()
        client_creds_file = config_dir / "client_credentials"
        session_tokens_file = config_dir / "session_tokens"

        # Check if credentials exist, if not setup
        if not client_creds_file.exists() or not session_tokens_file.exists():
            register_tokens(db)

        # Get credentials
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

        output_file = f"mlst_schemes_{db}.tab"
        processed_file = f"processed_dbs_{db}.tab"

        if not resume:
            clear_file(output_file)
            clear_file(processed_file)
        
        processed_dbs = load_processed_databases(processed_file)
        
        base_uri = BASE_API[db]
        resources = fetch_resources(base_uri, client_key, client_secret, 
                                  session_token, session_secret, verbose)
        
        if not resources:
            error("No resources found")
            sys.exit(1)

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

        success("Fetch completed successfully!")

    except Exception as e:
        error(f"An error occurred: {e}")
        info("Progress file kept for resume capability")
        if verbose:
            import traceback
            error(traceback.format_exc())
        sys.exit(1)