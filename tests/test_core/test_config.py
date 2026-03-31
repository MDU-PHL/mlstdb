import pytest
import stat
from unittest.mock import patch
from pathlib import Path
from mlstdb.core.config import check_dir


def test_check_dir_creates_missing_directory(tmp_path):
    """check_dir creates the directory when it does not exist."""
    new_dir = tmp_path / "new" / "nested"
    assert not new_dir.exists()
    check_dir(str(new_dir))
    assert new_dir.is_dir()


def test_check_dir_accepts_existing_writable_directory(tmp_path):
    """check_dir succeeds silently on an existing writable directory."""
    check_dir(str(tmp_path))  # should not raise


def test_check_dir_raises_on_unwritable_directory(tmp_path):
    """check_dir raises PermissionError when the directory is not writable."""
    ro_dir = tmp_path / "readonly"
    ro_dir.mkdir()
    ro_dir.chmod(stat.S_IRUSR | stat.S_IXUSR)  # r-x, no write
    try:
        with pytest.raises(PermissionError, match="Cannot write to directory"):
            check_dir(str(ro_dir))
    finally:
        # Restore permissions so tmp_path cleanup can remove the directory.
        ro_dir.chmod(stat.S_IRWXU)


def test_check_dir_raises_when_write_test_fails(tmp_path):
    """check_dir raises PermissionError when the actual write test raises OSError.

    This tests the NFS / network-filesystem scenario where os.access would
    return True but a real write still fails.
    """
    with patch("pathlib.Path.touch", side_effect=OSError("NFS write denied")):
        with pytest.raises(PermissionError, match="Cannot write to directory"):
            check_dir(str(tmp_path))


def test_check_dir_cleans_up_write_test_file(tmp_path):
    """check_dir leaves no .write_test artefact behind after a successful check."""
    check_dir(str(tmp_path))
    assert not (tmp_path / ".write_test").exists()
