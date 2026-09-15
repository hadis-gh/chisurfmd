#!/usr/bin/env bash

set -euo pipefail


# ==============================================================================
# INPUT
# ==============================================================================

ROOT_DIR="$1"
MD_EXE="$2"

OUTPUT_DIR="$ROOT_DIR/outputs/quickstart"
OUTPUT="$OUTPUT_DIR/quickstart_oriented_lj.bp"


# ==============================================================================
# CLEAN PREVIOUS OUTPUT
# ==============================================================================

rm -rf "$OUTPUT_DIR"


# ==============================================================================
# RUN QUICK START
# ==============================================================================

MD_EXE="$MD_EXE" \
bash "$ROOT_DIR/scripts/run_single.sh" \
     "$ROOT_DIR/scripts/configs/quickstart_oriented_lj.sh"


# ==============================================================================
# CHECK OUTPUT
# ==============================================================================

if [ ! -e "$OUTPUT" ]; then
    echo "FAILED: expected output was not created:"
    echo "  $OUTPUT"
    exit 1
fi

echo "Quick-start smoke test passed."