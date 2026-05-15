#!/usr/bin/env bash
# Sweep recycles x diffusion-steps for 5MHX (0.5occA / 0.5occB) across all
# four model wrappers. Per-model recycle lists skip each model's existing
# default-recycles baseline (which lives on a separate machine). Each
# (R, D) combo writes to its own --output-dir so runs do not collide
# (run_grid_search.py does not encode R or D in the per-run path).
set -euo pipefail

# ===================== User-configurable variables =====================
# Path to protein CSV consumed by run_grid_search.py --proteins.
PROTEINS_CSV="${PROTEINS_CSV:-5mhx_0.5occA_0.5occB.csv}"

# Base output dir; per-combo subdirs ${BASE_OUT}/recyc${R}_diff${D}/ are created.
BASE_OUT="${BASE_OUT:-grid_search_results_5mhx_recyc_diff}"

# Diffusion-steps axis (shared across all models).
DIFF_STEPS=(50 100 200 400)

# Per-model recycles lists. Default-recycles cell skipped for each model.
BOLTZ_RECYCLES=(1 10)        # Boltz default = 3
PROTENIX_RECYCLES=(1 3)      # Protenix default = 10 (model N_cycle)
RF3_RECYCLES=(1 3)           # RF3 default = 10

# Held-fixed pure_guidance config (matches existing native grid).
ENSEMBLE_SIZE=8
SCALER=pure_guidance
BOLTZ_GW=0.1
PROTENIX_GW=0.1
RF3_GW=0.01

# Extra flags passed to every invocation. Drop --dry-run when ready to run.
EXTRA_ARGS=()                # e.g. EXTRA_ARGS=(--dry-run) or (--force-all)
# =======================================================================

COMMON_ARGS=(
  --proteins "$PROTEINS_CSV"
  --scalers "$SCALER"
  --ensemble-sizes "$ENSEMBLE_SIZE"
  --augmentation
  --align-to-input
)

run_combo () {
  # $1 = pixi env, $2 = recycles, $3 = diffusion-steps, rest = extra model args
  local env="$1"; shift
  local r="$1"; shift
  local d="$1"; shift
  local out="${BASE_OUT}/recyc${r}_diff${d}"
  mkdir -p "$out"
  pixi run -e "$env" python run_grid_search.py \
    "${COMMON_ARGS[@]}" \
    --recycling-steps "$r" \
    --num-diffusion-steps "$d" \
    --output-dir "$out" \
    "$@" \
    "${EXTRA_ARGS[@]}"
}

# Boltz-2 (X-RAY DIFFRACTION + MD)
# for R in "${BOLTZ_RECYCLES[@]}"; do
#   for D in "${DIFF_STEPS[@]}"; do
#     run_combo boltz "$R" "$D" \
#       --model boltz2 --method "X-RAY DIFFRACTION" \
#       --gradient-weights "$BOLTZ_GW"
# 
#     run_combo boltz "$R" "$D" \
#       --model boltz2 --method MD \
#       --gradient-weights "$BOLTZ_GW"
#   done
# done

# Protenix
for R in "${PROTENIX_RECYCLES[@]}"; do
  for D in "${DIFF_STEPS[@]}"; do
    run_combo protenix "$R" "$D" \
      --model protenix \
      --gradient-weights "$PROTENIX_GW"
  done
done

# RF3
for R in "${RF3_RECYCLES[@]}"; do
  for D in "${DIFF_STEPS[@]}"; do
    run_combo rf3 "$R" "$D" \
      --model rf3 \
      --gradient-weights "$RF3_GW"
  done
done
