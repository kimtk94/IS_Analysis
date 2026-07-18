#!/usr/bin/env python3
"""Run the full restored-project QA audit."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEPS = [
    [sys.executable, "scripts/01_build_expected_manifest.py"],
    [sys.executable, "scripts/02_audit_raw_downloads.py"],
]


def main() -> int:
    for step in STEPS:
        print("Running:", " ".join(step), flush=True)
        result = subprocess.run(step, cwd=ROOT, check=False)
        if result.returncode:
            return result.returncode
    print("Full audit completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
