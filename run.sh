#!/usr/bin/env bash

# ------------------------------------------
# MCD UI Launcher
# Runs main.R using the newest available Rscript
# ------------------------------------------

set -u

# ------------------------------------------
# Configuration
# ------------------------------------------
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT_NAME="main.R"
MAIN_R="${SCRIPT_DIR}/${R_SCRIPT_NAME}"
RSCRIPT_EXE=""

# ------------------------------------------
# Check that the R script exists
# ------------------------------------------
if [[ ! -f "$MAIN_R" ]]; then
    echo "Error: Could not find ${R_SCRIPT_NAME} in:"
    echo "$SCRIPT_DIR"
    echo
    exit 1
fi

# ------------------------------------------
# Compare version strings
# returns 0 if $1 >= $2
# Handles versions like 4.3, 4.3.2, 4.4.0
# ------------------------------------------
version_ge() {
    local v1="$1"
    local v2="$2"

    local a1 b1 c1
    local a2 b2 c2

    IFS=. read -r a1 b1 c1 <<< "$v1"
    IFS=. read -r a2 b2 c2 <<< "$v2"

    b1="${b1:-0}"
    c1="${c1:-0}"
    b2="${b2:-0}"
    c2="${c2:-0}"

    a1="${a1//[^0-9]/}"
    b1="${b1//[^0-9]/}"
    c1="${c1//[^0-9]/}"
    a2="${a2//[^0-9]/}"
    b2="${b2//[^0-9]/}"
    c2="${c2//[^0-9]/}"

    a1="${a1:-0}"
    b1="${b1:-0}"
    c1="${c1:-0}"
    a2="${a2:-0}"
    b2="${b2:-0}"
    c2="${c2:-0}"

    if (( 10#$a1 > 10#$a2 )); then return 0; fi
    if (( 10#$a1 < 10#$a2 )); then return 1; fi

    if (( 10#$b1 > 10#$b2 )); then return 0; fi
    if (( 10#$b1 < 10#$b2 )); then return 1; fi

    if (( 10#$c1 >= 10#$c2 )); then return 0; fi

    return 1
}

# ------------------------------------------
# Get R version from an Rscript executable
# ------------------------------------------
get_r_version() {
    local rscript="$1"

    "$rscript" --vanilla -e 'cat(as.character(getRversion()))' 2>/dev/null
}

# ------------------------------------------
# Add candidate Rscript if executable
# ------------------------------------------
CANDIDATES=()

add_candidate() {
    local candidate="$1"

    if [[ -x "$candidate" ]]; then
        CANDIDATES+=("$candidate")
    fi
}

# ------------------------------------------
# Candidate locations
# ------------------------------------------

# Prefer PATH first, because this respects the user's active environment
if command -v Rscript >/dev/null 2>&1; then
    add_candidate "$(command -v Rscript)"
fi

# macOS R framework locations
add_candidate "/Library/Frameworks/R.framework/Resources/bin/Rscript"
add_candidate "/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript"

if [[ -d "/Library/Frameworks/R.framework/Versions" ]]; then
    for d in /Library/Frameworks/R.framework/Versions/*; do
        [[ -d "$d" ]] || continue
        add_candidate "$d/Resources/bin/Rscript"
        add_candidate "$d/bin/Rscript"
    done
fi

# Common Linux / HPC / custom locations
add_candidate "/usr/bin/Rscript"
add_candidate "/usr/local/bin/Rscript"
add_candidate "/opt/R/bin/Rscript"

if [[ -d "/opt/R" ]]; then
    for d in /opt/R/*; do
        [[ -d "$d" ]] || continue
        add_candidate "$d/bin/Rscript"
        add_candidate "$d/lib/R/bin/Rscript"
    done
fi

if [[ -d "/usr/local/lib/R" ]]; then
    add_candidate "/usr/local/lib/R/bin/Rscript"
fi

if [[ -d "/usr/lib/R" ]]; then
    add_candidate "/usr/lib/R/bin/Rscript"
fi

# ------------------------------------------
# Select newest working Rscript
# ------------------------------------------
BEST_VERSION=""

for candidate in "${CANDIDATES[@]}"; do
    [[ -x "$candidate" ]] || continue

    version="$(get_r_version "$candidate")"

    if [[ -z "$version" ]]; then
        continue
    fi

    if [[ -z "$RSCRIPT_EXE" ]] || version_ge "$version" "$BEST_VERSION"; then
        RSCRIPT_EXE="$candidate"
        BEST_VERSION="$version"
    fi
done

# ------------------------------------------
# If Rscript was not found, fail clearly
# ------------------------------------------
if [[ -z "$RSCRIPT_EXE" ]]; then
    echo "Error: Rscript was not found on this system."
    echo "Please install R and make sure Rscript is available."
    echo "See the README.md for setup instructions."
    echo
    exit 1
fi

# ------------------------------------------
# Run the R script
# ------------------------------------------
echo "=========================================="
echo "MCD UI Launcher"
echo "=========================================="
echo "Using Rscript:"
echo "$RSCRIPT_EXE"
echo
echo "Detected R version:"
echo "$BEST_VERSION"
echo
echo "Running script:"
echo "$MAIN_R"
echo "=========================================="
echo

"$RSCRIPT_EXE" --vanilla "$MAIN_R" "$@"
EXITCODE=$?

echo
if [[ "$EXITCODE" -ne 0 ]]; then
    echo "The script ended with an error."
else
    echo "The script finished successfully."
fi

exit "$EXITCODE"