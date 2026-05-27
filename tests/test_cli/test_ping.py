from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from mlstdb.cli.ping import ping


TEST_URL = "https://rest.pubmlst.org/db"
JSON_BODY = {"databases": [{"href": "https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef"}]}
PLAIN_BODY = "OK"


# ---------------------------------------------------------------------------
# Help text
# ---------------------------------------------------------------------------

def test_ping_help():
    runner = CliRunner()
    result = runner.invoke(ping, ["--help"])
    assert result.exit_code == 0
    assert "URL" in result.output
    assert "--no-auth" in result.output
    assert "--db" in result.output


# ---------------------------------------------------------------------------
# --no-auth path
# ---------------------------------------------------------------------------

def test_ping_no_auth_json_response():
    """--no-auth skips credentials; JSON response is pretty-printed."""
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.text = '{"databases": []}'
    mock_resp.json.return_value = {"databases": []}

    with patch("mlstdb.core.auth.requests.Session") as mock_session_cls:
        mock_session = MagicMock()
        mock_session.get.return_value = mock_resp
        mock_session_cls.return_value = mock_session

        runner = CliRunner()
        result = runner.invoke(ping, [TEST_URL, "--no-auth"])

    assert result.exit_code == 0
    assert "HTTP 200" in result.output
    assert "databases" in result.output


def test_ping_no_auth_plain_response():
    """--no-auth: plain-text response body is printed as-is."""
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.text = PLAIN_BODY
    mock_resp.json.side_effect = ValueError("not json")

    with patch("mlstdb.core.auth.requests.Session") as mock_session_cls:
        mock_session = MagicMock()
        mock_session.get.return_value = mock_resp
        mock_session_cls.return_value = mock_session

        runner = CliRunner()
        result = runner.invoke(ping, [TEST_URL, "--no-auth"])

    assert result.exit_code == 0
    assert PLAIN_BODY in result.output


# ---------------------------------------------------------------------------
# API key path
# ---------------------------------------------------------------------------

def test_ping_api_key_path():
    """When an API key is stored, it is used for authentication."""
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.text = "{}"
    mock_resp.json.return_value = {}

    with patch("mlstdb.core.auth.retrieve_api_key", return_value="testapikey123"), \
         patch("mlstdb.core.auth.requests.Session") as mock_session_cls:
        mock_session = MagicMock()
        mock_session.get.return_value = mock_resp
        mock_session_cls.return_value = mock_session

        runner = CliRunner()
        result = runner.invoke(ping, [TEST_URL, "--db", "pubmlst"])

    assert result.exit_code == 0
    assert "HTTP 200" in result.output
    # Confirm the X-API-Key header was set on the session
    mock_session.headers.update.assert_any_call(
        {"User-Agent": mock_session.headers.update.call_args_list[0][0][0]["User-Agent"],
         "X-API-Key": "testapikey123"}
    )


# ---------------------------------------------------------------------------
# OAuth path
# ---------------------------------------------------------------------------

def test_ping_oauth_path():
    """When no API key but OAuth tokens are present, OAuth is used."""
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.text = "{}"
    mock_resp.json.return_value = {}

    with patch("mlstdb.core.auth.retrieve_api_key", return_value=None), \
         patch("mlstdb.core.auth.get_client_credentials", return_value=("cid", "csecret")), \
         patch("mlstdb.core.auth.retrieve_session_token", return_value=("stoken", "ssecret")), \
         patch("mlstdb.core.auth.OAuth1Session") as mock_oauth_cls:
        mock_oauth = MagicMock()
        mock_oauth.get.return_value = mock_resp
        mock_oauth_cls.return_value = mock_oauth

        runner = CliRunner()
        result = runner.invoke(ping, [TEST_URL, "--db", "pubmlst"])

    assert result.exit_code == 0
    assert "HTTP 200" in result.output
    mock_oauth_cls.assert_called_once()


# ---------------------------------------------------------------------------
# Unauthenticated fallback (no creds stored)
# ---------------------------------------------------------------------------

def test_ping_unauthenticated_fallback():
    """When no creds are stored, ping falls back to unauthenticated with a warning."""
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.text = "{}"
    mock_resp.json.return_value = {}

    with patch("mlstdb.core.auth.retrieve_api_key", return_value=None), \
         patch("mlstdb.core.auth.get_client_credentials", side_effect=ValueError("no creds")), \
         patch("mlstdb.core.auth.requests.Session") as mock_session_cls:
        mock_session = MagicMock()
        mock_session.get.return_value = mock_resp
        mock_session_cls.return_value = mock_session

        runner = CliRunner()
        result = runner.invoke(ping, [TEST_URL, "--db", "pubmlst"])

    assert result.exit_code == 0
    assert "Warning" in result.output or "unauthenticated" in result.output.lower()
    assert "HTTP 200" in result.output


# ---------------------------------------------------------------------------
# 401 and 403 error handling
# ---------------------------------------------------------------------------

def test_ping_401_response():
    """401 responses print the status, body, and a short hint, then exit 0."""
    mock_resp = MagicMock()
    mock_resp.status_code = 401
    mock_resp.text = "Unauthorised"
    mock_resp.json.side_effect = ValueError("not json")

    with patch("mlstdb.core.auth.requests.Session") as mock_session_cls:
        mock_session = MagicMock()
        mock_session.get.return_value = mock_resp
        mock_session_cls.return_value = mock_session

        runner = CliRunner()
        result = runner.invoke(ping, [TEST_URL, "--no-auth"])

    assert result.exit_code == 0
    assert "HTTP 401" in result.output
    assert "mlstdb connect" in result.output


def test_ping_403_response():
    """403 responses print the status, body, and a short hint, then exit 0."""
    mock_resp = MagicMock()
    mock_resp.status_code = 403
    mock_resp.text = "Forbidden"
    mock_resp.json.side_effect = ValueError("not json")

    with patch("mlstdb.core.auth.requests.Session") as mock_session_cls:
        mock_session = MagicMock()
        mock_session.get.return_value = mock_resp
        mock_session_cls.return_value = mock_session

        runner = CliRunner()
        result = runner.invoke(ping, [TEST_URL, "--no-auth"])

    assert result.exit_code == 0
    assert "HTTP 403" in result.output
    assert "mlstdb connect" in result.output


# ---------------------------------------------------------------------------
# Missing URL argument
# ---------------------------------------------------------------------------

def test_ping_missing_url():
    """Omitting URL argument exits with error code 2."""
    runner = CliRunner()
    result = runner.invoke(ping, [])
    assert result.exit_code == 2
    assert "Missing argument" in result.output
