#!/usr/bin/env bash
# Copy the refined-patched ensemble CIFs (+ density inputs) for the cross-model
# side-chain guidance panel from the diffuse-shared mount into this folder with
# clean, PyMOL-safe names.
#
# Unlike the single-protein panels (8VPW, 5I09), this panel pairs ONE protein
# with ONE model each - the residue that best shows that model's side-chain
# guidance:
#   boltz2 X-ray : 3DT2 resi 11  (PHE)  gw0.0, gw0.1
#   boltz2 MD    : 2QCV resi 261 (LYS)  gw0.0, gw0.1
#   rf3          : 3L4P resi 18  (ASP)  gw0.0, gw0.01
#   protenix     : 3T94 resi 229 (THR)  gw0.0, gw0.1
#
# Variant copied : refined-patched.cif (the patched ensemble output)
# Occupancy      : 0.5occA_0.5occB
# Naming         : <PDB>_<model>_<gw>.cif  +  <PDB>_density_input.cif
#                  - dot-free stems => each filename is a valid PyMOL object name
#                  - gw decimal encoded with 'p'  (gw0p01 = guidance weight 0.01)
set -euo pipefail

ROOT="$HOME/astera/diffuse/diffuse-shared/raw/sampleworks"
DEST="$HOME/Downloads/sampleworks_paper/for_team/sidechain_panel"

[ -d "$ROOT/occ_sweep_results" ] || {
  echo "ERROR: source not found at $ROOT" >&2
  echo "       Is the diffuse-shared PVC mounted? (actl mount diffuse-shared)" >&2
  exit 1
}
mkdir -p "$DEST"

# Each entry: <path relative to ROOT>::<destination filename>
entries=(
  "occ_sweep_results/3DT2_0.5occA_0.5occB/boltz2_X-RAY_DIFFRACTION/pure_guidance/ens8_gw0.0/refined-patched.cif::3DT2_boltz2_xray_gw0p0.cif"
  "occ_sweep_results/3DT2_0.5occA_0.5occB/boltz2_X-RAY_DIFFRACTION/pure_guidance/ens8_gw0.1/refined-patched.cif::3DT2_boltz2_xray_gw0p1.cif"
  "initial_dataset_40_occ_sweeps/processed/3DT2/3DT2_single_001_density_input.cif::3DT2_density_input.cif"
  "occ_sweep_results/2QCV_0.5occA_0.5occB/boltz2_MD/pure_guidance/ens8_gw0.0/refined-patched.cif::2QCV_boltz2_md_gw0p0.cif"
  "occ_sweep_results/2QCV_0.5occA_0.5occB/boltz2_MD/pure_guidance/ens8_gw0.1/refined-patched.cif::2QCV_boltz2_md_gw0p1.cif"
  "initial_dataset_40_occ_sweeps/processed/2QCV/2QCV_single_001_density_input.cif::2QCV_density_input.cif"
  "occ_sweep_results_rf3_fixed/3L4P_0.5occA_0.5occB/rf3/pure_guidance/ens8_gw0.0/refined-patched.cif::3L4P_rf3_fixed_gw0p0.cif"
  "occ_sweep_results_rf3_fixed/3L4P_0.5occA_0.5occB/rf3/pure_guidance/ens8_gw0.01/refined-patched.cif::3L4P_rf3_fixed_gw0p01.cif"
  "initial_dataset_40_occ_sweeps/processed/3L4P/3L4P_single_001_density_input.cif::3L4P_density_input.cif"
  "occ_sweep_results/3T94_0.5occA_0.5occB/protenix/pure_guidance/ens8_gw0.0/refined-patched.cif::3T94_protenix_gw0p0.cif"
  "occ_sweep_results/3T94_0.5occA_0.5occB/protenix/pure_guidance/ens8_gw0.1/refined-patched.cif::3T94_protenix_gw0p1.cif"
  "initial_dataset_40_occ_sweeps/processed/3T94/3T94_single_001_density_input.cif::3T94_density_input.cif"
)

manifest="$DEST/MANIFEST.txt"
{
  echo "# Cross-model side-chain guidance panel - refined-patched ensembles + density inputs"
  echo "# One protein per model: 3DT2/boltz2-xray, 2QCV/boltz2-md, 3L4P/rf3, 3T94/protenix"
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
