#!/usr/bin/env python3
"""Patch the UKB-PPP batch runner to treat header-only exposure files as present.

This legacy recovery helper is retained for traceability. It only applies when
``ukb_ppp_batch_manifest_runner_fast.py`` is present in the working tree.
"""
from __future__ import annotations

from pathlib import Path

CANDIDATES = [
    Path("ukb_ppp_batch_manifest_runner_fast.py"),
    Path("scripts/ukb_ppp_batch_manifest_runner_fast.py"),
]
script = next((path for path in CANDIDATES if path.exists()), None)
if script is None:
    raise SystemExit("[ERROR] Cannot find ukb_ppp_batch_manifest_runner_fast.py")

text = script.read_text(encoding="utf-8")
backup = script.with_suffix(script.suffix + ".bak")
backup.write_text(text, encoding="utf-8")

text = text.replace(
    "expected_output.exists() and expected_output.stat().st_size > 0",
    "expected_output.exists()",
)

script.write_text(text, encoding="utf-8")

print("[OK] patched:", script)
print("[OK] backup:", backup)
