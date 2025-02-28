import configparser
import click
import sys
from pathlib import Path
from typing import Tuple
from rauth import OAuth1Service, OAuth1Session
from mlstdb.core.config import get_config_dir, BASE_API, BASE_WEB, DB_MAPPING
from mlstdb.utils import error, success, info

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
    
    # Message after registration
    click.secho("\n=== Registration Complete ===", fg="green", bold=True)
    click.echo("\nThe script will now fetch MLST scheme data, which may take a while.")
    click.echo("The process will download data from multiple databases and may take several minutes.")
    
    if not click.confirm("\nDo you want to continue with data fetching now?", default=True):
        info("You can run the script again later to fetch the data.")
        sys.exit(0)
    
    # If user wants to continue, return the tokens
    return token, secret

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