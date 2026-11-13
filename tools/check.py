"""Run from any directory: python /path/to/this/tools/check.py (no install)."""

import os
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
environment = os.environ.copy()
environment["PYTHONPATH"] = str(root / "src")
flags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
for arguments in (("-m", "unittest", "discover", "-s", "tests", "-v"),
                  ("-m", "fixedpoint_fir_lab", "experiment", "--output", "results/demo.json")):
    completed = subprocess.run([sys.executable, "-B", *arguments], cwd=root, env=environment,
                               capture_output=True, text=True, timeout=60, creationflags=flags)
    print(completed.stdout, end="", flush=True)
    print(completed.stderr, end="", file=sys.stderr, flush=True)
    completed.check_returncode()
