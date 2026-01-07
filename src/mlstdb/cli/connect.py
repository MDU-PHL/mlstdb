import click
import configparser
from mlstdb.core.auth import register_tokens, test_connection
from mlstdb.core.config import get_config_dir
from mlstdb.utils import error, success, info


@click.command()
@click.help_option('-h', '--help')
@click.option('--db', '-d', type=click.Choice(['pubmlst', 'pasteur']), 
              help='Database to use (pubmlst or pasteur)')
@click.option('--verbose', '-v', is_flag=True, 
              help='Enable verbose logging for debugging')
def connect(db, verbose):
    """Initial Database Registration and Setup
    
    Establishes connection with PubMLST or Pasteur databases by registering
    OAuth credentials and obtaining session tokens. This is required before
    using the update command.  
    """
    
    try:
        # If db is not provided, prompt for it
        if not db:
            db = click. prompt(
                "Which database would you like to connect to?",
                type=click.Choice(['pubmlst', 'pasteur']),
                default='pubmlst'
            )
        
        config_dir = get_config_dir()
        client_creds_file = config_dir / "client_credentials"
        session_tokens_file = config_dir / "session_tokens"
        
        # Check if already connected
        config = configparser.ConfigParser(interpolation=None)
        
        already_connected = False
        if client_creds_file.exists() and session_tokens_file.exists():
            config.read(client_creds_file)
            has_client_creds = config.has_section(db)
            
            config.read(session_tokens_file)
            has_session_tokens = config.has_section(db)
            
            if has_client_creds and has_session_tokens:
                already_connected = True
                click.secho(f"\n✓ Credentials found for {db}", fg="green")
                
                # Test the connection
                if test_connection(db, verbose=verbose):
                    success(f"\n✓ Connection to {db} is valid!")
                    info(f"Using existing credentials for {db}")
                    info("\nYou can now use 'mlstdb update' to update/download your database.")
                    return
                else:
                    error(f"\n✗ Connection test failed for {db}")
                    info("The credentials exist but the connection is not valid.")
                    
                    if not click. confirm(f"\nDo you want to re-register with {db}?", default=True):
                        error("\nmlstdb not connected. Please re-register when ready.") 
                        raise SystemExit(1)
                    # If user says yes, continue to registration below
        
        # Register tokens
        if verbose:
            info(f"Starting registration process for {db}...")
        
        register_tokens(db)
        
        # Test the newly registered connection
        info("\nVerifying new connection...")
        if test_connection(db, verbose=verbose):
            success(f"\n✓ Successfully connected to {db}!")
            info("\nNext steps:")
            info("  1. Use 'mlstdb update' to update/download MLST schemes")
            info("  2. Or use 'mlstdb fetch' for advanced schema exploration")
        else:
            error(f"\n✗ Connection test failed after registration")
            info("Please check your credentials and try again.")
            raise SystemExit(1)
        
    except KeyboardInterrupt:
        error("\n\nConnection cancelled by user")
        raise SystemExit(1)
    except Exception as e: 
        error(f"Connection failed: {e}")
        if verbose:
            import traceback
            error(traceback.format_exc())
        raise SystemExit(1)