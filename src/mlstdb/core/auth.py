import configparser
import click
import os
import sys
import requests
from pathlib import Path
from typing import Tuple, Optional
from rauth import OAuth1Service, OAuth1Session
from mlstdb.core.config import get_config_dir, BASE_API, BASE_WEB, DB_MAPPING
from mlstdb.utils import error, success, info
from mlstdb.__about__ import __version__


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
    
    fd = os.open(file_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as configfile:
        config.write(configfile)
    success(f"\nClient credentials saved to {file_path}")
    return client_id, client_secret


def setup_api_key(site: str) -> str:
    """Prompt for and save a personal API key (BIGSdb ≥ v1.53.0)."""
    config = configparser.ConfigParser(interpolation=None)
    file_path = get_config_dir() / "api_keys"
    if file_path.exists():
        config.read(file_path)

    info("\nPlease enter your personal API key (from your BIGSdb profile page):")
    api_key = click.prompt("API Key", type=str).strip()
    while not api_key:
        error("API key cannot be empty")
        api_key = click.prompt("API Key", type=str).strip()

    config[site] = {"api_key": api_key}
    fd = os.open(file_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as configfile:
        config.write(configfile)
    success(f"\nAPI key saved to {file_path}")
    return api_key


def retrieve_api_key(site: str) -> Optional[str]:
    """Return the stored personal API key for *site*, or None if not found."""
    config = configparser.ConfigParser(interpolation=None)
    file_path = get_config_dir() / "api_keys"
    if file_path.is_file():
        config.read(file_path)
        if config.has_section(site):
            return config[site].get("api_key")
    return None


def _parse_error_message(response) -> str:
    """Safely extract an error message from a failed API response."""
    try:
        data = response.json()
        if isinstance(data, dict):
            return data.get("message", f"HTTP {response.status_code}")
        return f"HTTP {response.status_code}"
    except Exception:
        text = response.text[:200] if response.text else "(empty response)"
        return f"HTTP {response.status_code} — {text}"


def register_tokens(db: str):
    """Setup authentication tokens by registering with the service."""
    info(f"\nNo tokens found for {db}. Starting registration process...")
    
    # Setup client credentials
    client_id, client_secret = setup_client_credentials(db)
    
    # Initialize OAuth service
    service = OAuth1Service(
        name="MLSTdb downloader",
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
        headers={"User-Agent": f"mlstdb/{__version__}"}
    )
    if r.status_code != 200:
        msg = _parse_error_message(r)
        if "timestamp" in msg.lower() or "600 seconds" in msg.lower():
            error(
                f"Failed to get request token: {msg}\n"
                "  This is a clock synchronisation issue.\n"
                "  On WSL/Linux, try:  sudo hwclock -s\n"
                "  On WSL2, try:       sudo ntpdate time.windows.com"
            )
        else:
            error(f"Failed to get request token: {msg}")
        sys.exit(1)
    
    request_token = r.json()["oauth_token"]
    request_secret = r.json()["oauth_token_secret"]
    success("Temporary token received")
    
    # Get access token
    click.secho("\nAuthorisation Required", fg="yellow", bold=True)
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
        headers={"User-Agent": f"mlstdb/{__version__}"},
    )
    
    if r.status_code != 200:
        error(f"Failed to get access token: {_parse_error_message(r)}")
        sys.exit(1)
        
    access_token = r. json()["oauth_token"]
    access_secret = r.json()["oauth_token_secret"]
    
    # Save access token
    config = configparser.ConfigParser(interpolation=None)
    file_path = get_config_dir() / "access_tokens"
    if file_path.exists():
        config.read(file_path)
    config[db] = {"token": access_token, "secret": access_secret}
    fd = os.open(file_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as configfile:
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
    session.headers.update({"User-Agent": f"mlstdb/{__version__}"})

    r = session.get(url, params={})

    if r.status_code != 200:
        error(f"Failed to get session token: {_parse_error_message(r)}")
        sys.exit(1)
        
    token = r.json()["oauth_token"]
    secret = r. json()["oauth_token_secret"]
    
    # Save session token
    config = configparser. ConfigParser(interpolation=None)
    file_path = get_config_dir() / "session_tokens"
    if file_path.exists():
        config.read(file_path)
    config[db] = {"token": token, "secret": secret}
    fd = os.open(file_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as configfile:
        config.write(configfile)
    
    success(f"\nSession token saved to {file_path}")
    
    # Message after registration
    click.secho("\n=== Registration Complete ===", fg="green", bold=True)
    
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
    for file_name in ["client_credentials", "session_tokens", "access_tokens", "api_keys"]:
        file_path = config_dir / file_name
        if file_path.exists():
            config = configparser.ConfigParser(interpolation=None)
            config.read(file_path)
            if db in config:
                config. remove_section(db)
                with open(file_path, 'w') as f:
                    config.write(f)
                success(f"Removed {db} credentials from {file_name}")

def retrieve_session_token(key_name: str) -> Tuple[str, str]:
    """Get OAuth session token from config file."""
    config = configparser.ConfigParser(interpolation=None)
    file_path = get_config_dir() / "session_tokens"
    
    if file_path. is_file():
        config.read(file_path)
        if config.has_section(key_name):
            return (config[key_name]["token"], 
                   config[key_name]["secret"])
    
    return None, None

def test_connection(db: str, verbose: bool = False, api_key: Optional[str] = None) -> bool:
    """Test if the connection to the database is valid.

    Args:
        db: Database name ('pubmlst' or 'pasteur')
        verbose: If True, display JSON payload from test URI
        api_key: Optional personal API key (BIGSdb ≥ v1.53.0)

    Returns:
        True if connection is valid, False otherwise
    """
    try:
        # Test URL - using the database info endpoint
        test_url = f"{BASE_API[db]}/db/{DB_MAPPING[db]}/schemes"

        info(f"\nTesting connection to {db}...")
        info(f"Using test database:  {DB_MAPPING[db]}")

        if api_key:
            # API key auth path (BIGSdb ≥ v1.53.0)
            session = requests.Session()
            session.headers.update({
                "User-Agent": f"mlstdb/{__version__}",
                "X-API-Key": api_key,
            })
            if verbose:
                info(f"\nRequesting:  {test_url}")
            response = session.get(test_url)
            if verbose:
                info(f"\nResponse status code: {response.status_code}")
                if response.status_code == 200:
                    try:
                        import json as json_module
                        info("\nJSON payload received:")
                        click.echo(json_module.dumps(response.json(), indent=2))
                    except Exception as e:
                        error(f"Could not parse JSON response: {e}")
            if response.status_code == 200:
                return True
            elif response.status_code == 401:
                error("\nAPI key rejected (401). The key may be revoked or invalid.")
                info("Run 'mlstdb connect --api-key' to save a new API key.")
                return False
            else:
                error(f"\nConnection test failed with status code: {response.status_code}")
                return False

        # OAuth path
        info("\nPlease ensure you are registered to this database.")

        # Get client credentials
        client_id, client_secret = get_client_credentials(db)

        # Get session tokens
        session_token, session_secret = retrieve_session_token(db)

        if not session_token or not session_secret:
            return False

        # Create OAuth session
        session = OAuth1Session(
            consumer_key=client_id,
            consumer_secret=client_secret,
            access_token=session_token,
            access_token_secret=session_secret,
        )
        session.headers.update({"User-Agent": f"mlstdb/{__version__}"})

        # Make test request
        if verbose:
            info(f"\nRequesting:  {test_url}")

        response = session.get(test_url, params={})
        
        if verbose:
            info(f"\nResponse status code: {response.status_code}")
            if response.status_code == 200:
                try:
                    json_payload = response.json()
                    info("\nJSON payload received:")
                    import json as json_module
                    click.echo(json_module.dumps(json_payload, indent=2))
                except Exception as e:
                    error(f"Could not parse JSON response: {e}")
        
        if response.status_code == 200:
            return True
        elif response.status_code == 401:
            # Try to refresh the session token
            info("\nAttempting to refresh session token...")
            new_tokens = refresh_session_token(db, client_id, client_secret, verbose)
            if new_tokens:
                info("Session token refreshed. Testing connection again...")
                # Retry the test with new token
                session_token, session_secret = new_tokens
                session = OAuth1Session(
                    consumer_key=client_id,
                    consumer_secret=client_secret,
                    access_token=session_token,
                    access_token_secret=session_secret,
                )
                session.headers. update({"User-Agent": f"mlstdb/{__version__}"})
                response = session.get(test_url, params={})
                
                if response. status_code == 200:
                    success("Connection successful after token refresh!")
                    return True
                else:
                    error(f"\nAuthentication failed after token refresh - status code: {response.status_code}")
                    return False
            else:
                error("\nAuthentication failed - session token may be invalid or expired")
                return False
        elif response.status_code == 403:
            error("\nAccess denied - you may not have permission to access this database")
            info(f"Please ensure you are registered to the '{DB_MAPPING[db]}' database")
            return False
        else:
            error(f"\nConnection test failed with status code: {response. status_code}")
            return False
            
    except ValueError as e:
        error(f"\n{e}")
        return False
    except Exception as e:
        error(f"\nConnection test failed: {e}")
        if verbose:
            import traceback
            error(traceback.format_exc())
        return False

def refresh_session_token(db:  str, client_key: str, client_secret: str, verbose: bool = False) -> Optional[Tuple[str, str]]: 
    """Refresh session token using existing access tokens. 
    
    Args:
        db: Database name ('pubmlst' or 'pasteur')
        client_key: OAuth client ID
        client_secret: OAuth client secret
        verbose: If True, display detailed information
        
    Returns: 
        Tuple of (new_token, new_secret) if successful, None otherwise
    """
    try:
        info("Invalid session token. Requesting new one...")
        
        # Get access tokens
        config = configparser.ConfigParser(interpolation=None)
        access_tokens_file = get_config_dir() / "access_tokens"
        
        if not access_tokens_file.exists():
            error("Access tokens file not found. Please re-register.")
            return None
        
        config.read(access_tokens_file)
        if not config.has_section(db):
            error(f"No access tokens found for {db}. Please re-register.")
            return None
            
        access_token = config[db]["token"]
        access_secret = config[db]["secret"]
        
        # Get new session token
        url_session = f"{BASE_API[db]}/db/{DB_MAPPING[db]}/oauth/get_session_token"
        session_request = OAuth1Session(
            client_key,
            client_secret,
            access_token=access_token,
            access_token_secret=access_secret,
        )
        session_request.headers.update({"User-Agent": f"mlstdb/{__version__}"})
        
        if verbose:
            info(f"Requesting new session token from:  {url_session}")
        
        r = session_request.get(url_session, params={})
        
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
                config. write(configfile)
            
            if verbose:
                success("New session token obtained and saved")
            else:
                success("Session token refreshed successfully")
            
            return new_token, new_secret
        else:
            error(f"Failed to refresh session token: {r.status_code}")
            if verbose and r.text:
                error(f"Response: {r.text}")
            return None
            
    except Exception as e:
        error(f"Failed to refresh session token: {e}")
        if verbose:
            import traceback
            error(traceback.format_exc())


def ping_url(
    url: str,
    db: Optional[str] = None,
    verbose: bool = False,
    no_auth: bool = False,
) -> tuple:
    """Make a single GET request to *url* and return the result.

    Auth priority (unless *no_auth* is True):
    1. API key stored for *db* (``X-API-Key`` header).
    2. OAuth session token stored for *db*.
    3. Unauthenticated fallback with a warning.

    Args:
        url: The API URL to probe.
        db: Database name (``'pubmlst'`` or ``'pasteur'``).  Required for
            authenticated paths; ignored when *no_auth* is ``True``.
        verbose: Print extra request detail.
        no_auth: Skip all authentication and use a plain session.

    Returns:
        Tuple of ``(status_code: int, body: str, json_payload: dict | None)``.
        *json_payload* is ``None`` when the response body is not valid JSON.
    """
    import json as _json

    session = None
    auth_mode = "unauthenticated"

    if no_auth:
        session = requests.Session()
        session.headers.update({"User-Agent": f"mlstdb/{__version__}"})
        auth_mode = "no-auth"
    else:
        # 1. Try API key
        api_key = retrieve_api_key(db) if db else None
        if api_key:
            session = requests.Session()
            session.headers.update({
                "User-Agent": f"mlstdb/{__version__}",
                "X-API-Key": api_key,
            })
            auth_mode = "api-key"
        else:
            # 2. Try OAuth session token
            try:
                client_id, client_secret = get_client_credentials(db)
                session_token, session_secret = retrieve_session_token(db)
                if session_token and session_secret:
                    session = OAuth1Session(
                        consumer_key=client_id,
                        consumer_secret=client_secret,
                        access_token=session_token,
                        access_token_secret=session_secret,
                    )
                    session.headers.update({"User-Agent": f"mlstdb/{__version__}"})
                    auth_mode = "oauth"
            except ValueError:
                pass

            # 3. Unauthenticated fallback
            if session is None:
                click.secho(
                    f"\nWarning: No credentials found for '{db}', trying unauthenticated...",
                    fg="yellow",
                )
                session = requests.Session()
                session.headers.update({"User-Agent": f"mlstdb/{__version__}"})
                auth_mode = "unauthenticated"

    if verbose:
        info(f"\nAuth mode:   {auth_mode}")
        info(f"Requesting:  {url}")

    response = session.get(url)

    if verbose:
        info(f"Status:      {response.status_code}")

    body = response.text
    json_payload = None
    try:
        json_payload = response.json()
    except Exception:
        pass

    return response.status_code, body, json_payload