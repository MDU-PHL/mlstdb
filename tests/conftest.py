import pytest
from pathlib import Path

@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / "test_data"

@pytest.fixture
def scheme_uris_file(test_data_dir):
    return test_data_dir / "scheme_uris.tab"