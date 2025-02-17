import pytest
from click.testing import CliRunner
from mlstdb.cli.fetch import fetch

def test_fetch_command_help():
    runner = CliRunner()
    result = runner.invoke(fetch, ['--help'])
    assert result.exit_code == 0
    assert 'Database to use' in result.output

def test_fetch_with_invalid_db():
    runner = CliRunner()
    result = runner.invoke(fetch, ['-d', 'invalid'])
    assert result.exit_code == 2
    assert 'Invalid value for' in result.output