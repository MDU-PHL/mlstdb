import click
from pathlib import Path
import importlib.resources
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import configparser
import sys
import os
from mlstdb.core.auth import register_tokens, setup_client_credentials, remove_db_credentials
from mlstdb.core.download import (fetch_resources, get_matching_schemes, create_session,
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
@click.option('--no-auth', is_flag=True,
              help='Use unauthenticated access (faster, works for public APIs)')
@click.option('--threads', '-t', default=1, show_default=True,
              help='Number of parallel threads for fetching schemes, max: 4 to avoid too many simultaneous requests')
@click.option('--verbose', '-v', is_flag=True, 
              help='Enable verbose logging for debugging')
def fetch(db, exclude, match, scheme_uris, filter, resume, no_auth, threads, verbose):
    """[ADVANCED] BIGSdb Scheme Explorer and Fetcher
    
    ⚠️  ADVANCED USE ONLY - For scheme exploration and custom workflows. 
    
    For standard usage, use the recommended workflow:
      1. mlstdb connect --db pubmlst
      2. mlstdb connect --db pasteur
      3. mlstdb update
    
    The fetch subcommand allows exploration of all available schemes from
    BIGSdb databases with custom filtering and matching options.
    It handles both authentication and scheme discovery.
    """
    
    # Show deprecation notice
    click.secho("\n" + "="*70, fg="yellow")
    click.secho("  ADVANCED COMMAND - Not required for standard usage", fg="yellow", bold=True)
    click.secho("="*70, fg="yellow")
    click.echo("\nFor most users, please use the recommended workflow:")
    click.secho("  1. mlstdb connect --db pubmlst/pasteur", fg="cyan")
    click.secho("  2. mlstdb update", fg="cyan")
    click.echo("\nThis 'fetch' command is for:")
    click.echo("  • Exploring all available schemes")
    click.echo("  • Custom filtering of schemes")
    click.secho("\n" + "="*70 + "\n", fg="yellow")
    
    if not click.confirm("Do you want to continue with advanced fetch?", default=False):
        info("Cancelled. Use 'mlstdb connect' and 'mlstdb update' instead.")
        sys.exit(0)
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
        
        client_key = None
        client_secret = None
        session_token = None
        session_secret = None

        if no_auth:
            info("Using unauthenticated access (no OAuth)")
            http_session = create_session(no_auth=True)
        else:
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

            # Create a single reusable OAuth session
            http_session = create_session(client_key, client_secret, 
                                          session_token, session_secret)

        output_file = f"mlst_schemes_{db}.tab"
        processed_file = f"processed_dbs_{db}.tab"

        if not resume:
            clear_file(output_file)
            clear_file(processed_file)
        
        processed_dbs = load_processed_databases(processed_file)
        
        base_uri = BASE_API[db]
        resources = fetch_resources(base_uri, client_key, client_secret, 
                                  session_token, session_secret, verbose,
                                  session=http_session)
        
        if not resources:
            error("No resources found")
            sys.exit(1)

        # Collect all databases to process
        all_databases = []
        for resource in resources:
            if 'databases' in resource:
                for database in resource['databases']:
                    if database['description'] in processed_dbs:
                        if verbose:
                            info(f"Skipping already processed database: {database['description']}")
                        continue
                    all_databases.append(database)

        total_dbs = len(all_databases) + len(processed_dbs)
        success_count = len(processed_dbs)
        auth_skipped = []
        all_scheme_lines = []

        # Process databases in parallel
        info(f"Processing {len(all_databases)} databases using {threads} threads...")
        
        def _process_db(database):
            """Worker function for parallel database processing."""
            return get_matching_schemes(
                database, match, exclude,
                client_key, client_secret,
                session_token, session_secret,
                verbose=verbose,
                session=http_session,
            )

        with ThreadPoolExecutor(max_workers=threads) as executor:
            future_to_db = {
                executor.submit(_process_db, database): database
                for database in all_databases
            }
            for future in tqdm(as_completed(future_to_db), 
                              total=len(all_databases),
                              desc="Processing databases"):
                database = future_to_db[future]
                try:
                    result = future.result()
                    if result.get('auth_skipped'):
                        auth_skipped.append(result['auth_skipped'])
                    if result.get('matching_schemes'):
                        all_scheme_lines.extend(result['matching_schemes'])
                    success_count += 1
                    # Track progress for resume
                    save_processed_database(processed_file, database['description'])
                except Exception as e:
                    error(f"Error processing {database['description']}: {e}")
                    continue

        # Write all matching schemes to the output file at once
        if all_scheme_lines:
            with open(output_file, 'a') as f:
                f.writelines(all_scheme_lines)

        if auth_skipped:
            click.secho(
                "\nThe following databases were skipped due to authentication issues:",
                fg="yellow"
            )
            for db_name in auth_skipped:
                click.secho(f"  - {db_name}", fg="yellow")
            click.secho(
                "\nPlease ensure you are registered in the database and have the necessary "
                "permissions to access these schemes. Run 'mlstdb connect' to set up credentials.",
                fg="yellow"
            )

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

        success("Fetch completed successfully! View the results in " + output_file + "\n" + "Use `mlstdb update` to download the required MLST datasets.")

    except Exception as e:
        error(f"An error occurred: {e}")
        info("Progress file kept for resume capability")
        if verbose:
            import traceback
            error(traceback.format_exc())
        sys.exit(1)