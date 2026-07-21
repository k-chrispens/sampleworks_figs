#!/usr/bin/env bash
# Copy the 8VPW / 0.5occA_0.5occB refined-patched ensemble CIFs (+ density input)
# from the diffuse-shared mount into this folder with clean, PyMOL-safe names.
#
# Variant copied : refined-patched.cif (the patched ensemble output)
# Occupancy      : 0.5occA_0.5occB
# Models         : the 4 MAIN models only (boltz2 X-ray, boltz2 MD, protenix, rf3)
# Guidance       : baseline vs guided (two states per model)
#                  - boltz2/protenix : gw0.0, gw0.1
#                  - rf3             : gw0.0, gw0.01
# Naming         : 8VPW_<model>_<gw>.cif
#                  - dot-free stems => each filename is a valid PyMOL object name
#                  - gw decimal encoded with 'p'  (gw0p01 = guidance weight 0.01)
set -euo pipefail

ROOT="$HOME/astera/diffuse/diffuse-shared/raw/sampleworks"
DEST="$HOME/Downloads/sampleworks_paper/for_team/8VPW_results"
OCC="8VPW_0.5occA_0.5occB"

[ -d "$ROOT/occ_sweep_results" ] || {
  echo "ERROR: source not found at $ROOT" >&2
  echo "       Is the diffuse-shared PVC mounted? (actl mount diffuse-shared)" >&2
  exit 1
}
mkdir -p "$DEST"

# Each entry: <path relative to ROOT>::<destination filename>
entries=(
  "occ_sweep_results/$OCC/boltz2_X-RAY_DIFFRACTION/pure_guidance/ens8_gw0.0/refined-patched.cif::8VPW_boltz2_xray_gw0p0.cif"
  "occ_sweep_results/$OCC/boltz2_X-RAY_DIFFRACTION/pure_guidance/ens8_gw0.1/refined-patched.cif::8VPW_boltz2_xray_gw0p1.cif"
  "occ_sweep_results/$OCC/boltz2_MD/pure_guidance/ens8_gw0.0/refined-patched.cif::8VPW_boltz2_md_gw0p0.cif"
  "occ_sweep_results/$OCC/boltz2_MD/pure_guidance/ens8_gw0.1/refined-patched.cif::8VPW_boltz2_md_gw0p1.cif"
  "occ_sweep_results/$OCC/protenix/pure_guidance/ens8_gw0.0/refined-patched.cif::8VPW_protenix_gw0p0.cif"
  "occ_sweep_results/$OCC/protenix/pure_guidance/ens8_gw0.1/refined-patched.cif::8VPW_protenix_gw0p1.cif"
  "occ_sweep_results_rf3_fixed/$OCC/rf3/pure_guidance/ens8_gw0.0/refined-patched.cif::8VPW_rf3_fixed_gw0p0.cif"
  "occ_sweep_results_rf3_fixed/$OCC/rf3/pure_guidance/ens8_gw0.01/refined-patched.cif::8VPW_rf3_fixed_gw0p01.cif"
  "initial_dataset_40_occ_sweeps/processed/8VPW/8VPW_single_001_density_input.cif::8VPW_density_input.cif"
)

manifest="$DEST/MANIFEST.txt"
{
  echo "# 8VPW 0.5occA_0.5occB - refined-patched ensembles + density input (4 main models)"
  echo "# dest_filename  <-  source (relative to $ROOT)"
} > "$manifest"

ok=0; fail=0
for e in "${entries[@]}"; do
  src="${e%%::*}"; dst="${e##*::}"
  if [ -f "$ROOT/$src" ]; then
    cp -f "$ROOT/$src" "$DEST/$dst"
    printf '%-32s  <-  %s\n' "$dst" "$src" >> "$manifest"
    printf 'OK    %s\n' "$dst"
    ok=$((ok + 1))
  else
    printf 'MISS  %s  (source missing: %s)\n' "$dst" "$src" >&2
    fail=$((fail + 1))
  fi
done

echo
echo "Copied $ok file(s) into $DEST  (failures: $fail)"
echo "Manifest written to $manifest"
