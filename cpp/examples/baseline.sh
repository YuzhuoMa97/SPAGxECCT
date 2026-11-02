#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/.."
make test
python3 tests/run_integration.py
