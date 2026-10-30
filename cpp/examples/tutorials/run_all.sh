#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/../.."
bash examples/tutorials/01_spagxe.sh
bash examples/tutorials/02_spagxe_plus.sh
bash examples/tutorials/03_spagxemix.sh
bash examples/tutorials/04_spagxemix_local.sh
bash examples/tutorials/05_global_local_cct.sh
bash examples/tutorials/06_local_packing.sh
