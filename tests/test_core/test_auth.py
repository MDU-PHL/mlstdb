import pytest
from unittest.mock import patch
from mlstdb.core.auth import register_tokens, setup_client_credentials, setup_api_key, retrieve_api_key
from pathlib import Path


@pytest.fixture(autouse=True)
def clear_mlstdb_env_vars(monkeypatch):
    """Ensure tests in this module run isolated from host environment variables by default."""
    for env in [
        "MLSTDB_PUBMLST_API_KEY",
        "MLSTDB_PASTEUR_API_KEY",
        "MLSTDB_PUBMLST_CLIENT_ID",
        "MLSTDB_PUBMLST_CLIENT_SECRET",
        "MLSTDB_PASTEUR_CLIENT_ID",
        "MLSTDB_PASTEUR_CLIENT_SECRET",
    ]:
        monkeypatch.delenv(env, raising=False)

def test_setup_client_credentials(tmp_path):
    config_dir = tmp_path / ".bigsdb_tokens"
    config_dir.mkdir()
    
    # Mock both get_config_dir and click.prompt
    with patch('mlstdb.core.auth.get_config_dir', return_value=config_dir), \
         patch('click.prompt', side_effect=["123456" * 4, "1234567" * 6]):  # 24 and 42 chars
        result = setup_client_credentials("pubmlst")
        assert result is not None
        assert (config_dir / "client_credentials").exists()
        
        # Verify the content of the credentials file
        with open(config_dir / "client_credentials", 'r') as f:
            content = f.read()
            assert "123456" * 4 in content  # 24-char client ID
            assert "1234567" * 6 in content  # 42-char secret


def test_setup_client_credentials_uses_env_vars(monkeypatch):
    monkeypatch.setenv("MLSTDB_PASTEUR_CLIENT_ID", "123456789012345678901234")
    monkeypatch.setenv("MLSTDB_PASTEUR_CLIENT_SECRET", "123456789012345678901234567890123456789012")
    client_id, client_secret = setup_client_credentials("pasteur")
    assert client_id == "123456789012345678901234"
    assert client_secret == "123456789012345678901234567890123456789012"


def test_setup_api_key(tmp_path):
    config_dir = tmp_path / ".config" / "mlstdb"
    config_dir.mkdir(parents=True)

    with patch('mlstdb.core.auth.get_config_dir', return_value=config_dir), \
         patch('click.prompt', return_value='my-test-api-key-abc123'):
        result = setup_api_key('pubmlst')

    assert result == 'my-test-api-key-abc123'
    api_keys_file = config_dir / 'api_keys'
    assert api_keys_file.exists()
    import stat
    assert oct(stat.S_IMODE(api_keys_file.stat().st_mode)) == '0o600'
    with open(api_keys_file) as f:
        content = f.read()
    assert 'my-test-api-key-abc123' in content
    assert '[pubmlst]' in content


def test_retrieve_api_key_found(tmp_path):
    config_dir = tmp_path / ".config" / "mlstdb"
    config_dir.mkdir(parents=True)

    api_keys_file = config_dir / 'api_keys'
    api_keys_file.write_text('[pubmlst]\napi_key = abc-xyz-123\n')

    with patch('mlstdb.core.auth.get_config_dir', return_value=config_dir):
        key = retrieve_api_key('pubmlst')

    assert key == 'abc-xyz-123'


def test_retrieve_api_key_missing(tmp_path, monkeypatch):
    monkeypatch.delenv("MLSTDB_PUBMLST_API_KEY", raising=False)
    config_dir = tmp_path / ".config" / "mlstdb"
    config_dir.mkdir(parents=True)

    with patch('mlstdb.core.auth.get_config_dir', return_value=config_dir):
        key = retrieve_api_key('pubmlst')

    assert key is None


def test_retrieve_api_key_from_env(monkeypatch):
    monkeypatch.setenv("MLSTDB_PUBMLST_API_KEY", "env-api-key-123")
    assert retrieve_api_key("pubmlst") == "env-api-key-123"


def test_retrieve_api_key_env_precedence_over_file(tmp_path, monkeypatch):
    config_dir = tmp_path / ".config" / "mlstdb"
    config_dir.mkdir(parents=True)
    api_keys_file = config_dir / 'api_keys'
    api_keys_file.write_text('[pubmlst]\napi_key = file-key-456\n')

    monkeypatch.setenv("MLSTDB_PUBMLST_API_KEY", "env-key-789")

    with patch('mlstdb.core.auth.get_config_dir', return_value=config_dir):
        assert retrieve_api_key('pubmlst') == "env-key-789"


def test_get_client_credentials_from_env(monkeypatch):
    from mlstdb.core.auth import get_client_credentials
    monkeypatch.setenv("MLSTDB_PUBMLST_CLIENT_ID", "env_client_id_1234567890")
    monkeypatch.setenv("MLSTDB_PUBMLST_CLIENT_SECRET", "env_client_secret_12345678901234567890")

    client_id, client_secret = get_client_credentials("pubmlst")
    assert client_id == "env_client_id_1234567890"
    assert client_secret == "env_client_secret_12345678901234567890"


def test_get_client_credentials_env_precedence_over_file(tmp_path, monkeypatch):
    from mlstdb.core.auth import get_client_credentials
    config_dir = tmp_path / ".config" / "mlstdb"
    config_dir.mkdir(parents=True)
    creds_file = config_dir / 'client_credentials'
    creds_file.write_text('[pubmlst]\nclient_id = file_id\nclient_secret = file_secret\n')

    monkeypatch.setenv("MLSTDB_PUBMLST_CLIENT_ID", "env_id")
    monkeypatch.setenv("MLSTDB_PUBMLST_CLIENT_SECRET", "env_secret")

    with patch('mlstdb.core.auth.get_config_dir', return_value=config_dir):
        cid, csec = get_client_credentials("pubmlst")
        assert cid == "env_id"
        assert csec == "env_secret"


def test_get_client_credentials_missing(tmp_path, monkeypatch):
    from mlstdb.core.auth import get_client_credentials
    monkeypatch.delenv("MLSTDB_PUBMLST_CLIENT_ID", raising=False)
    monkeypatch.delenv("MLSTDB_PUBMLST_CLIENT_SECRET", raising=False)
    config_dir = tmp_path / ".config" / "mlstdb"
    config_dir.mkdir(parents=True)

    with patch('mlstdb.core.auth.get_config_dir', return_value=config_dir):
        with pytest.raises(ValueError, match="Client credentials not found for pubmlst"):
            get_client_credentials("pubmlst")