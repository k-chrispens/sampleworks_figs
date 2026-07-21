#!/usr/bin/env bash
# Copy the inputs for the 3-structure "examples" panel (examples_panel.pml) into
# this folder with clean, flat, PyMOL-safe names so the script is self-contained.
#
# Structures : 6B8X, 5MHX, 5SOP   (the unified set from baseline/view.pml)
# Per structure, 4 MAIN models, each baseline (gw0.0) vs guided (optimal gw):
#     boltz2 X-ray  gw0.0 / gw0.1   (good_darkblue)
#     boltz2 MD     gw0.0 / gw0.1   (good_teal)
#     protenix      gw0.0 / gw0.1   (good_pink)
#     rf3           gw0.0 / gw0.01  (good_navaho)
# Variant copied : refined-patched.cif (the patched ens8 ensemble; 8 states)
# Occupancy      : 0.5occA_0.5occB
# Plus per structure: ground-truth density input (gray reference, keeps alt A/B)
#                     and the experimental 1.00 A map (ccp4, for isomesh blobs).
#
# Two sources:
#   ROOT   - diffuse-shared PVC mount: refined-patched CIFs + density inputs
#   MAPSRC - local baseline run dir  : the .ccp4 maps (mount has no maps)
#
# Naming: dot-free stems => each filename is a valid PyMOL object name and the
#         filename matches the object name used in examples_panel.pml 1:1.
#         gw decimal encoded with 'p' (gw0p1 = guidance weight 0.1).
set -euo pipefail

ROOT="$HOME/astera/diffuse/diffuse-shared/raw/sampleworks"
MAPSRC="$HOME/Downloads/initial_dataset_40/outputs/initial_dataset_40/baseline"
DEST="$HOME/Downloads/sampleworks_paper/for_team/examples_panels"

[ -d "$ROOT/occ_sweep_results" ] || {
  echo "ERROR: source not found at $ROOT" >&2
  echo "       Is the diffuse-shared PVC mounted? (actl mount diffuse-shared)" >&2
  exit 1
}
mkdir -p "$DEST"

STRUCTS=(6B8X 5MHX 5SOP)

manifest="$DEST/MANIFEST.txt"
{
  echo "# examples_panels - inputs for examples_panel.pml"
  echo "# 3 structures (6B8X, 5MHX, 5SOP) x 4 main models, baseline gw0.0 vs guided optimal gw"
  echo "# refined-patched ens8 CIFs + density input from: $ROOT"
  echo "# experimental .ccp4 maps from: $MAPSRC"
  echo "# dest_filename  <-  source"
} > "$manifest"

ok=0; fail=0
copy() {  # copy <abs-src> <dest-basename> <manifest-src-label>
  local src="$1" dst="$2" label="$3"
  if [ -f "$src" ]; then
    cp -f "$src" "$DEST/$dst"
    printf '%-30s  <-  %s\n' "$dst" "$label" >> "$manifest"
    printf 'OK    %s\n' "$dst"
    ok=$((ok + 1))
  else
    printf 'MISS  %s  (source missing: %s)\n' "$dst" "$src" >&2
    fail=$((fail + 1))
  fi
}

for S in "${STRUCTS[@]}"; do
  OCC="${S}_0.5occA_0.5occB"

  # ground-truth density input (gray reference; keeps alt A/B) - from mount
  copy "$ROOT/initial_dataset_40_occ_sweeps/processed/$S/${S}_single_001_density_input.cif" \
       "${S}_gt.cif" "initial_dataset_40_occ_sweeps/processed/$S/${S}_single_001_density_input.cif"

  # experimental 1.00 A map (isomesh) - local only
  copy "$MAPSRC/${S}_0.5occA_0.5occB_1.00A.ccp4" \
       "${S}_map.ccp4" "(local baseline) ${S}_0.5occA_0.5occB_1.00A.ccp4"

  # boltz2 X-ray : gw0.0 / gw0.1
  copy "$ROOT/occ_sweep_results/$OCC/boltz2_X-RAY_DIFFRACTION/pure_guidance/ens8_gw0.0/refined-patched.cif" \
       "${S}_boltz2XRAY_gw0p0.cif" "occ_sweep_results/$OCC/boltz2_X-RAY_DIFFRACTION/.../ens8_gw0.0/refined-patched.cif"
  copy "$ROOT/occ_sweep_results/$OCC/boltz2_X-RAY_DIFFRACTION/pure_guidance/ens8_gw0.1/refined-patched.cif" \
       "${S}_boltz2XRAY_gw0p1.cif" "occ_sweep_results/$OCC/boltz2_X-RAY_DIFFRACTION/.../ens8_gw0.1/refined-patched.cif"

  # boltz2 MD : gw0.0 / gw0.1
  copy "$ROOT/occ_sweep_results/$OCC/boltz2_MD/pure_guidance/ens8_gw0.0/refined-patched.cif" \
       "${S}_boltz2MD_gw0p0.cif" "occ_sweep_results/$OCC/boltz2_MD/.../ens8_gw0.0/refined-patched.cif"
  copy "$ROOT/occ_sweep_results/$OCC/boltz2_MD/pure_guidance/ens8_gw0.1/refined-patched.cif" \
       "${S}_boltz2MD_gw0p1.cif" "occ_sweep_results/$OCC/boltz2_MD/.../ens8_gw0.1/refined-patched.cif"

  # protenix : gw0.0 / gw0.1
  copy "$ROOT/occ_sweep_results/$OCC/protenix/pure_guidance/ens8_gw0.0/refined-patched.cif" \
       "${S}_protenix_gw0p0.cif" "occ_sweep_results/$OCC/protenix/.../ens8_gw0.0/refined-patched.cif"
  copy "$ROOT/occ_sweep_results/$OCC/protenix/pure_guidance/ens8_gw0.1/refined-patched.cif" \
       "${S}_protenix_gw0p1.cif" "occ_sweep_results/$OCC/protenix/.../ens8_gw0.1/refined-patched.cif"

  # rf3 (rf3_fixed) : gw0.0 / gw0.01
  copy "$ROOT/occ_sweep_results_rf3_fixed/$OCC/rf3/pure_guidance/ens8_gw0.0/refined-patched.cif" \
       "${S}_rf3_gw0p0.cif" "occ_sweep_results_rf3_fixed/$OCC/rf3/.../ens8_gw0.0/refined-patched.cif"
  copy "$ROOT/occ_sweep_results_rf3_fixed/$OCC/rf3/pure_guidance/ens8_gw0.01/refined-patched.cif" \
       "${S}_rf3_gw0p01.cif" "occ_sweep_results_rf3_fixed/$OCC/rf3/.../ens8_gw0.01/refined-patched.cif"
done

echo
echo "Copied $ok file(s) into $DEST  (failures: $fail)"
echo "Manifest written to $manifest"
