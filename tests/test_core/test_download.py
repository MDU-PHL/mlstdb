import pytest
from mlstdb.core.download import create_session


def test_create_session_api_key():
    session = create_session(api_key='test-key-abc')
    assert session.headers.get('X-API-Key') == 'test-key-abc'
    assert 'mlstdb/' in session.headers.get('User-Agent', '')


def test_create_session_no_auth():
    session = create_session(no_auth=True)
    assert 'X-API-Key' not in session.headers
    assert 'mlstdb/' in session.headers.get('User-Agent', '')


def test_create_session_api_key_takes_precedence_over_no_auth():
    """api_key should win when both api_key and no_auth are set."""
    session = create_session(no_auth=True, api_key='should-win')
    assert session.headers.get('X-API-Key') == 'should-win'
