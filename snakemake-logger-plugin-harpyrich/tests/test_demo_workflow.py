"""
Test the demo workflow using the rich logger plugin.

This test runs the demo Snakefile to ensure the rich logger plugin works correctly.
"""

import subprocess
import tempfile
from pathlib import Path
import pytest


def test_demo_workflow():
    """Test that a sample workflow runs successfully with the rich logger plugin."""

    project_root = Path(__file__).parent.parent
    test_snakefile = project_root / "tests" / "fixtures" / "Snakefile"

    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "demo_output"
        output_dir.mkdir()

        cmd = [
            "snakemake",
            "-s",
            str(test_snakefile),
            "-d",
            str(output_dir),
            "--show-failed-logs",
            "--printshellcmds",
            "--logger",
            "rich",
            "--cores",
            "1",
            "output.txt",
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,  # 5 minute timeout
            )

            assert result.returncode == 0, (
                f"Snakemake failed with return code {result.returncode}.\n"
                f"STDOUT: {result.stdout}\n"
                f"STDERR: {result.stderr}"
            )

        except subprocess.TimeoutExpired:
            pytest.fail("Snakemake command timed out after 5 minutes")
