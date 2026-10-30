#!/usr/bin/env bash
# Works as a bundle-root launcher and from cpp/examples in a source checkout.
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
if [ -f "$script_dir/cpp/examples/tutorials/run_all.sh" ]; then
    cpp_dir="$script_dir/cpp"
elif [ -f "$script_dir/tutorials/run_all.sh" ]; then
    cpp_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)
else
    echo "Cannot locate cpp/examples/tutorials/run_all.sh; extract the complete bundle." >&2
    exit 1
fi
if [ -f "$cpp_dir/SPAGxE.exe" ]; then
    exe="$cpp_dir/SPAGxE.exe"
elif [ -f "$cpp_dir/SPAGxE" ]; then
    exe="$cpp_dir/SPAGxE"
elif [ -n "${SPAGXE_BIN:-}" ]; then
    exe="$SPAGXE_BIN"
elif [ -f "$cpp_dir/build/bin/SPAGxE.exe" ]; then
    exe="$cpp_dir/build/bin/SPAGxE.exe"
else
    exe="$cpp_dir/build/bin/SPAGxE"
fi
if [ ! -f "$exe" ]; then
    echo "SPAGxE executable is missing: $exe" >&2
    exit 1
fi
# A bundle always selects its own executable; resolve source-checkout overrides
# before tutorials change their working directory.
exe_dir=$(CDPATH= cd -- "$(dirname -- "$exe")" && pwd)
export SPAGXE_BIN="$exe_dir/$(basename -- "$exe")"
exec bash "$cpp_dir/examples/tutorials/run_all.sh"
