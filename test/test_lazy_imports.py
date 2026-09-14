"""Regression tests for imports that should remain lazy."""

import subprocess
import sys


def assert_modules_not_loaded(import_statement, module_names):
    """Run an import in a clean interpreter and inspect loaded modules."""
    code = "\n".join(
        [
            "import sys",
            import_statement,
            f"module_names = {module_names!r}",
            "loaded = [name for name in module_names if name in sys.modules]",
            "if loaded:",
            "    raise RuntimeError(f'Unexpected imports: {loaded}')",
        ]
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        check=False,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_frame_does_not_import_io():
    assert_modules_not_loaded("import pycbc.frame", ("pycbc.io",))


def test_gracedb_does_not_import_plotting_modules():
    assert_modules_not_loaded(
        "import pycbc.io.gracedb",
        ("pycbc.results", "matplotlib.pyplot"),
    )
