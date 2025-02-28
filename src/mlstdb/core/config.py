from pathlib import Path
import os
from typing import Dict

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

def get_config_dir() -> Path:
    """Create and return the configuration directory."""
    config_dir = Path.home() / ".config" / "mlstdb"
    if not config_dir.exists():
        config_dir.mkdir(parents=True, mode=0o700)
    return config_dir

def check_dir(directory: str) -> None:
    """Ensure the directory exists and is writable."""
    path = Path(directory)
    if not path.exists():
        path.mkdir(parents=True)
    if not (path.is_dir() and os.access(directory, os.W_OK)):
        raise PermissionError(f"Cannot write to directory: {directory}")