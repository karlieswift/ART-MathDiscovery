#!/usr/bin/env python3
from __future__ import annotations

import ast
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCAN_ROOT = Path("\\\\?\\" + str(ROOT)) if os.name == "nt" else ROOT
FORBIDDEN_TOP_LEVEL = {"docs", "paper_source_data", "provenance"}
FORBIDDEN_SUFFIXES = {".log", ".tex", ".aux", ".toc", ".synctex", ".fls", ".fdb_latexmk"}
FORBIDDEN_DIRECTORY_NAMES = {"log", "logs", "output", "outputs", "raw", "validation", "deliverable", "provenance"}


def main() -> int:
    errors: list[str] = []

    for name in sorted(FORBIDDEN_TOP_LEVEL):
        if (ROOT / name).exists():
            errors.append(f"forbidden top-level path: {name}")

    files = [path for path in SCAN_ROOT.rglob("*") if path.is_file() and ".git" not in path.parts]
    for path in files:
        relative = path.relative_to(SCAN_ROOT)
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            errors.append(f"forbidden artifact: {relative}")
        if any(part.lower() in FORBIDDEN_DIRECTORY_NAMES for part in relative.parts[:-1]):
            errors.append(f"forbidden artifact directory: {relative}")

    python_files = [path for path in files if path.suffix == ".py"]
    for path in python_files:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (OSError, SyntaxError, UnicodeError) as exc:
            errors.append(f"invalid Python source {path.relative_to(SCAN_ROOT)}: {exc}")

    print(f"[CODE-ONLY] files={len(files)} python={len(python_files)} errors={len(errors)}")
    for error in errors[:100]:
        print(f"[CODE-ONLY] {error}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
