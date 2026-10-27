#!/usr/bin/env python3
"""Run a byte-preserved original checker with optimization explicitly rejected."""

if not __debug__:
    raise SystemExit("Optimized Python is unsupported; run without -O/-OO.")

import argparse
from pathlib import Path
import runpy
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("calculation", choices=("focus", "diagrams"))
    args = parser.parse_args()
    names = {"focus": "check-focus.py", "diagrams": "check-local-diagrams.py"}
    target = (Path(__file__).resolve().parents[1] / "inputs/original-source/scripts"
              / names[args.calculation])
    sys.argv = [str(target)]
    runpy.run_path(str(target), run_name="__main__")


if __name__ == "__main__":
    main()
