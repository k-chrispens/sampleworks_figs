#!/usr/bin/env bash
# Pull CSV results from the remote occ_sweep_results* folders, tagging each
# copied file with the model/folder it came from.
#
# For each folder under $REMOTE_BASE we list depth-1 *.csv files, keep those
# whose filename matches $REGEX, and rsync each one down renamed as
#   <stem>_<model-label>.csv
# so e.g. occ_sweep_results_protenix_mini/rscc_results.csv arrives locally as
#   rscc_results_protenix_mini.csv
#
# The model label is derived from the folder name (the part after
# "occ_sweep_results"); the bare base folder is labelled "boltz_protenix".
#
# Existing local files are never clobbered silently: each one prompts
# [y/N] (default skip) before being overwritten. Set FORCE=1 to overwrite
# without asking.
#
# Usage:
#   REMOTE=me@host ./pull_remote_results.sh '<regex>'
#   REMOTE=me@host DEST_DIR=initial40results ./pull_remote_results.sh 'rscc|lddt'
#   REMOTE=me@host DRY_RUN=1 ./pull_remote_results.sh '.*'   # preview, no transfer
#   REMOTE=me@host FORCE=1   ./pull_remote_results.sh '.*'   # overwrite, no prompt
set -euo pipefail

# ===================== User-configurable variables =====================
# ssh target for the machine holding the data. Override on the command line:
#   REMOTE=me@host ./pull_remote_results.sh ...
REMOTE="${REMOTE:-user@host}"

# Base directory on the remote containing the occ_sweep_results* folders.
REMOTE_BASE="${REMOTE_BASE:-/mnt/diffuse-shared/raw/sampleworks}"

# Where pulled CSVs land locally.
DEST_DIR="${DEST_DIR:-.}"

# Regex matched (via grep -E) against each candidate file's basename.
# Defaults to "." = every depth-1 *.csv. Pass your own as the first argument.
REGEX="${1:-.}"

# Folders to sweep (depth 1, *.csv only).
FOLDERS=(
  occ_sweep_results
  occ_sweep_results_protenix_mini
  occ_sweep_results_protenix_tiny
  occ_sweep_results_rf3_chiral_off
  occ_sweep_results_rf3_fixed
)

# rsync flags for real transfers. DRY_RUN=1 is handled in the loop below by
# skipping the transfer entirely rather than via `rsync --dry-run`: macOS's
# openrsync DOES honour --dry-run, but still prints "Transfer starting" and
# the filename, which reads exactly like a real copy. Skipping rsync gives
# unambiguous "[DRY RUN]" output instead.
RSYNC_OPTS=(-avP)
# =======================================================================

if [ "$REMOTE" = "user@host" ]; then
  echo "ERROR: set REMOTE to your ssh target, e.g. REMOTE=me@host $0 '$REGEX'" >&2
  exit 1
fi

# Derive the model label appended to files from a given folder.
#   occ_sweep_results                  -> boltz_protenix
#   occ_sweep_results_protenix_mini    -> protenix_mini
#   occ_sweep_results_rf3_chiral_off   -> rf3_chiral_off
model_label () {
  local suffix="${1#occ_sweep_results}"   # strip common prefix
  suffix="${suffix#_}"                     # strip leading underscore, if any
  echo "${suffix:-boltz_protenix}"        # bare base folder -> boltz_protenix
}

mkdir -p "$DEST_DIR"
echo "Remote : $REMOTE:$REMOTE_BASE"
echo "Regex  : $REGEX"
echo "Dest   : $DEST_DIR"
echo

pulled=0
skipped=0
for folder in "${FOLDERS[@]}"; do
  label="$(model_label "$folder")"

  # List depth-1 *.csv files on the remote, filter by regex on the basename.
  matches="$(
    ssh "$REMOTE" "find '$REMOTE_BASE/$folder' -maxdepth 1 -type f -name '*.csv'" 2>/dev/null \
      | while IFS= read -r path; do
          bn="${path##*/}"
          printf '%s\n' "$bn" | grep -Eq "$REGEX" && printf '%s\n' "$path"
        done || true
  )"

  if [ -z "$matches" ]; then
    echo "[$folder -> $label] no matching csv files"
    continue
  fi

  while IFS= read -r path; do
    [ -z "$path" ] && continue
    bn="${path##*/}"          # e.g. rscc_results.csv
    stem="${bn%.csv}"         # e.g. rscc_results
    dest="$DEST_DIR/${stem}_${label}.csv"

    # Dry run: report what would happen and transfer nothing.
    if [ -n "${DRY_RUN:-}" ]; then
      if [ -e "$dest" ]; then
        echo "[DRY RUN] $bn  ->  $dest  (would overwrite existing)"
      else
        echo "[DRY RUN] $bn  ->  $dest"
      fi
      pulled=$((pulled + 1))
      continue
    fi

    # Warn before clobbering an existing local file (FORCE=1 overwrites
    # without asking). Read the reply from /dev/tty, not the loop's stdin,
    # which is the here-string feeding `path`.
    if [ -e "$dest" ] && [ -z "${FORCE:-}" ]; then
      if [ -r /dev/tty ]; then
        printf 'WARNING: %s already exists. Overwrite? [y/N] ' "$dest" >&2
        read -r reply < /dev/tty || reply=""
        case "$reply" in
          [Yy]|[Yy][Ee][Ss]) ;;
          *) echo "  skipped." >&2; skipped=$((skipped + 1)); continue ;;
        esac
      else
        echo "WARNING: $dest exists; skipping (set FORCE=1 to overwrite)." >&2
        skipped=$((skipped + 1))
        continue
      fi
    fi

    echo "[$folder -> $label] ${bn}  ->  ${dest}"
    rsync "${RSYNC_OPTS[@]}" "$REMOTE:$path" "$dest"
    pulled=$((pulled + 1))
  done <<< "$matches"
done

echo
if [ -n "${DRY_RUN:-}" ]; then
  echo "Dry run: ${pulled} file(s) would be pulled from ${#FOLDERS[@]} folders. Nothing was transferred."
else
  echo "Done. ${pulled} file(s) pulled, ${skipped} skipped from ${#FOLDERS[@]} folders."
fi
