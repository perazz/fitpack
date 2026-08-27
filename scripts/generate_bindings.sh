#!/bin/bash
# Regenerate the C/C++ bindings from fitpack_bindings.yaml.
#
# Usage (from anywhere — the script locates the repository from its own path):
#   ./scripts/generate_bindings.sh [--dry-run] [--verbose]
#
# The generator is the fortran-arrays binding generator, which this repository does NOT
# vendor: point FITPACK_BINDINGS_GENERATOR at your checkout, or keep fortran-arrays as a
# sibling of the repository and the default below finds it.
#
#   FITPACK_BINDINGS_GENERATOR=/path/to/fortran-arrays/tools/bindings/generate-bindings \
#       ./scripts/generate_bindings.sh
#
# Requirements: python3 with fparser, jinja2, click and pyyaml (see the generator's
# requirements.txt). A .venv next to the generator is activated automatically if present.
#
# Everything this writes is a GENERATED file — never hand-edit one. To change the C/C++
# surface, edit the Fortran source, fitpack_bindings.yaml, or a snippet under extra_methods/,
# then re-run this script. The hand-written parts of the interface are preserved here:
#
#   src/capi/fitpack_core_c.f90  the 25 fp_*_c core procedural bindings + 2 helpers
#   include/fitpack_core_c.h     their pure-C header
#   include/fpPoint.hpp          the fixed-dimension point type

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG="$REPO_ROOT/fitpack_bindings.yaml"
SRC_OUT="$REPO_ROOT/src/capi"
INC_OUT="$REPO_ROOT/include"

GENERATOR="${FITPACK_BINDINGS_GENERATOR:-$REPO_ROOT/../fortran-arrays/tools/bindings/generate-bindings}"

if [ ! -f "$GENERATOR" ]; then
    echo "error: binding generator not found at '$GENERATOR'." >&2
    echo "       Set FITPACK_BINDINGS_GENERATOR to your fortran-arrays checkout, e.g." >&2
    echo "       FITPACK_BINDINGS_GENERATOR=/path/to/fortran-arrays/tools/bindings/generate-bindings" >&2
    exit 1
fi

GENERATOR_DIR="$(cd "$(dirname "$GENERATOR")" && pwd)"
if [ -f "$GENERATOR_DIR/.venv/bin/activate" ]; then
    # shellcheck disable=SC1091
    source "$GENERATOR_DIR/.venv/bin/activate"
fi
export PYTHONPATH="$GENERATOR_DIR${PYTHONPATH:+:$PYTHONPATH}"

# Hand-written files that live in the same directories as the generated ones.
KEEP_F90="fitpack_core_c.f90"
KEEP_H="fitpack_core_c.h"
KEEP_HPP="fpPoint.hpp"

echo "=== FITPACK binding generation ==="
echo "Generator: $GENERATOR"
echo "Config:    $CONFIG"
echo "Fortran:   $SRC_OUT"
echo "Headers:   $INC_OUT"
echo ""

mkdir -p "$SRC_OUT"

# Drop the previous output so a type removed from the config cannot leave a stale file behind.
find "$SRC_OUT" -maxdepth 1 -name '*_c.f90'       ! -name "$KEEP_F90" -delete
find "$SRC_OUT" -maxdepth 1 -name '*_c_types.f90'                     -delete
rm -f "$SRC_OUT/fitpack_fx_status.f90"
find "$INC_OUT" -maxdepth 1 -name '*.h'           ! -name "$KEEP_H"   -delete
find "$INC_OUT" -maxdepth 1 -name '*.hpp'         ! -name "$KEEP_HPP" -delete

"$GENERATOR" --config "$CONFIG" --output-dir "$SRC_OUT" --standalone "$@"

# The generator also emits shared-library export-control scripts (.def/.ver/.exports). This
# repository ships a static archive built by fpm, which uses none of them.
rm -rf "$INC_OUT/export_control"

echo ""
echo "=== Generated ==="
echo "src/capi/:"
ls "$SRC_OUT"/*_c.f90 "$SRC_OUT"/*_c_types.f90 2>/dev/null | sed 's|.*/|  |' || echo "  (none)"
echo "include/:"
ls "$INC_OUT"/*.h "$INC_OUT"/*.hpp 2>/dev/null | sed 's|.*/|  |' || echo "  (none)"
echo ""
echo "Done. Run 'fpm test' to compile and check."
