import pytest
from unittest.mock import patch

from mlstdb.__about__ import __version__
from mlstdb.core.download import fetch_json

def test_user_agent_header():
    captured_headers = {}

    def mock_get(self, url, **kwargs):
        # Capture headers sent in the request
        captured_headers.update(self.headers)
        class MockResponse:
            status_code = 200
            def json(self):
                return {}
            def raise_for_status(self):
                pass
        return MockResponse()

    # Patch OAuth1Session.get to mock network calls
    with patch("rauth.OAuth1Session.get", new=mock_get):
        # You must supply dummy credentials and a URL
        fetch_json(
            url="https://pubmlst.org/api/db",
            client_key="dummy_key",
            client_secret="dummy_secret",
            session_token="dummy_token",
            session_secret="dummy_session_secret",
            verbose=False
        )
    # Check that the User-Agent header is as expected    
    expected_ua = f"mlstdb/{__version__}"
    assert captured_headers['User-Agent'] == expected_ua