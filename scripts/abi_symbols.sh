#!/bin/bash
# Extract the C ABI of the Fortran binding layer, and diff it against another git revision.
#
# Usage (from anywhere — the script locates the repository from its own path):
#   ./scripts/abi_symbols.sh                 # list the current symbols, one `name(args)` per line
#   ./scripts/abi_symbols.sh --at <ref>      # list the symbols of a git revision instead
#   ./scripts/abi_symbols.sh --diff <ref>    # classify the current ABI against <ref> (e.g. main, 1.0.0)
#
# A symbol is any procedure carrying `bind(C, name='...')` under src/. The argument list comes
# from the Fortran declaration, so the diff separates a pure rename from a signature change.
# Types declared `bind(C)` are storage layout, not symbols, and are not listed. Both sides of a
# diff scan src/ recursively, so a revision predating the src/capi/ layout compares cleanly.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

AWK_EXTRACT='
function flush_logical(line,   name, rest, depth, i, ch, args, start) {
    if (line !~ /[Bb][Ii][Nn][Dd][ \t]*\([ \t]*[Cc][ \t]*,/) return
    if (line !~ /[Nn][Aa][Mm][Ee][ \t]*=/) return

    # Exported C name, from the bind(C, name=...) clause.
    rest = line
    sub(/^.*[Bb][Ii][Nn][Dd][ \t]*\([ \t]*[Cc][ \t]*,[ \t]*[Nn][Aa][Mm][Ee][ \t]*=[ \t]*/, "", rest)
    if (rest ~ /^"/)      { sub(/^"/, "", rest); sub(/".*$/, "", rest) }
    else if (rest ~ /^'"'"'/) { sub(/^'"'"'/, "", rest); sub(/'"'"'.*$/, "", rest) }
    else return
    name = rest

    # Dummy-argument list, from the parenthesis that follows the procedure name.
    start = match(line, /[Ss][Uu][Bb][Rr][Oo][Uu][Tt][Ii][Nn][Ee][ \t]+[A-Za-z_][A-Za-z_0-9]*[ \t]*\(/)
    if (start == 0) start = match(line, /[Ff][Uu][Nn][Cc][Tt][Ii][Oo][Nn][ \t]+[A-Za-z_][A-Za-z_0-9]*[ \t]*\(/)
    if (start == 0) { print name "()"; return }

    i = start + RLENGTH        # first character after the opening parenthesis
    depth = 1
    args = ""
    while (i <= length(line) && depth > 0) {
        ch = substr(line, i, 1)
        if (ch == "(") depth++
        else if (ch == ")") { depth--; if (depth == 0) break }
        if (depth > 0) args = args ch
        i++
    }
    gsub(/[ \t]+/, "", args)
    print tolower(name) "(" tolower(args) ")"
}
{
    line = $0
    # Drop a trailing comment, but only where no quote could be hiding the "!".
    if (line !~ /["'"'"']/) sub(/!.*$/, "", line)
    sub(/[ \t]+$/, "", line)
    if (cont != "") { sub(/^[ \t]*&/, "", line); line = cont line; cont = "" }
    if (line ~ /&$/) { sub(/&$/, "", line); cont = line; next }
    flush_logical(line)
}
END { if (cont != "") flush_logical(cont) }
'

# Symbols of the working tree.
current_symbols() {
    # shellcheck disable=SC2046
    awk "$AWK_EXTRACT" $(find "$REPO_ROOT/src" -name '*_c.f90' | LC_ALL=C sort) | LC_ALL=C sort -u
}

# Symbols of a git revision.
ref_symbols() {
    local ref="$1" f
    for f in $(git -C "$REPO_ROOT" ls-tree -r --name-only "$ref" -- src | grep '_c\.f90$'); do
        git -C "$REPO_ROOT" show "$ref:$f" | awk "$AWK_EXTRACT"
    done | LC_ALL=C sort -u
}

if [ "${1:-}" = "--at" ]; then
    ref_symbols "${2:?usage: $0 --at <git-ref>}"
    exit 0
fi

if [ "${1:-}" != "--diff" ]; then
    current_symbols
    exit 0
fi

REF="${2:?usage: $0 --diff <git-ref>}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

ref_symbols "$REF"  > "$TMP/old"
current_symbols     > "$TMP/new"

# name -> signature, for the rename/signature classification.
cut -d'(' -f1 "$TMP/old" | LC_ALL=C sort -u > "$TMP/old.names"
cut -d'(' -f1 "$TMP/new" | LC_ALL=C sort -u > "$TMP/new.names"

comm -12 "$TMP/old.names" "$TMP/new.names" > "$TMP/kept.names"
comm -23 "$TMP/old.names" "$TMP/new.names" > "$TMP/removed.names"
comm -13 "$TMP/old.names" "$TMP/new.names" > "$TMP/added.names"

: > "$TMP/identical"
: > "$TMP/changed"
while read -r n; do
    o="$(grep -m1 "^$n(" "$TMP/old")"
    w="$(grep -m1 "^$n(" "$TMP/new")"
    if [ "$o" = "$w" ]; then echo "$o" >> "$TMP/identical"
    else printf '%s\n    old: %s\n    new: %s\n' "$n" "${o#*(}" "${w#*(}" >> "$TMP/changed"
    fi
done < "$TMP/kept.names"

n_old=$(wc -l < "$TMP/old")
n_new=$(wc -l < "$TMP/new")
n_same=$(wc -l < "$TMP/identical")
n_chg=$(wc -l < "$TMP/kept.names")
n_chg=$((n_chg - n_same))
n_rm=$(wc -l < "$TMP/removed.names")
n_add=$(wc -l < "$TMP/added.names")

echo "# ABI diff: $REF -> working tree"
echo
echo "| class | count |"
echo "|-------|-------|"
printf '| symbols at %s | %s |\n' "$REF" "$n_old"
printf '| symbols now | %s |\n' "$n_new"
printf '| identical (name + signature) | %s |\n' "$n_same"
printf '| same name, changed signature | %s |\n' "$n_chg"
printf '| removed | %s |\n' "$n_rm"
printf '| new | %s |\n' "$n_add"
echo
echo "## Same name, changed signature"
echo
cat "$TMP/changed"
echo
echo "## Removed"
echo
cat "$TMP/removed.names"
echo
echo "## New"
echo
cat "$TMP/added.names"
