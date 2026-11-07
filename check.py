#!/usr/bin/env python3
"""Run offline checks and report passed, failed and skipped cases."""
import sys
import unittest
from pathlib import Path

if __name__ == "__main__":
    root=Path(__file__).resolve().parent
    loader=unittest.TestLoader()
    suite=loader.discover(str(root),pattern="test_configuration.py")
    suite.addTests(loader.discover(str(root/"tests"),pattern="test_*.py",top_level_dir=str(root/"tests")))
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    print(f"Passed: {result.testsRun-len(result.failures)-len(result.errors)-len(result.skipped)}; "
          f"failed: {len(result.failures)+len(result.errors)}; skipped: {len(result.skipped)}")
    sys.exit(0 if result.wasSuccessful() else 1)
