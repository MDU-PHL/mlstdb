import pytest
from unittest.mock import patch
from mlstdb.core.auth import register_tokens, setup_client_credentials
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