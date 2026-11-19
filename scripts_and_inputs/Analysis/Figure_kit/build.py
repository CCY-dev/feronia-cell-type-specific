#!/usr/bin/env python
"""One-command build for the FER figure kit. Run from this folder:  python build.py   (see README.md)."""
import os
import sys
from pathlib import Path

os.environ.setdefault('SOURCE_DATE_EPOCH', '1767225600')   # 2026-01-01: fixed timestamps in outputs
os.environ.setdefault('MPLBACKEND', 'Agg')
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ferfigs.cli import main  # noqa: E402

if __name__ == '__main__':
    sys.exit(main())
