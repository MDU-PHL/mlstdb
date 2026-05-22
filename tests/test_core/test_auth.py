import pytest
from unittest.mock import patch
from mlstdb.core.auth import register_tokens, setup_client_credentials, setup_api_key, retrieve_api_key
from pathlib import Path

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


def test_retrieve_api_key_missing(tmp_path):
    config_dir = tmp_path / ".config" / "mlstdb"
    config_dir.mkdir(parents=True)

    with patch('mlstdb.core.auth.get_config_dir', return_value=config_dir):
        key = retrieve_api_key('pubmlst')

    assert key is None