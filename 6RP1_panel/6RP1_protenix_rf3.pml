# 6RP1: density input + protenix + RF3 predictions
# Render protenix and rf3 separately by toggling show/hide before each ray.

# Custom color for rf3 (matches MODEL_PALETTE in style.py)
set_color rf3_gold, [1.0, 0.843, 0.0]

# Density input (reference)
load ./6RP1_single_001_density_input.cif, density_input
color gray70, density_input
util.cnc density_input
show cartoon, density_input

# Protenix predictions
load ./protenix_6RP1_0.0.cif, protenix_baseline
align protenix_baseline, density_input
color good_pink, protenix_baseline
util.cnc protenix_baseline
show cartoon, protenix_baseline

load ./protenix_6RP1_0.1.cif, protenix_gw_0_1
align protenix_gw_0_1, density_input
color good_pink, protenix_gw_0_1
util.cnc protenix_gw_0_1
show cartoon, protenix_gw_0_1

# RF3 predictions
load ./rf3_6RP1_0.0.cif, rf3_baseline
align rf3_baseline, density_input
color rf3_gold, rf3_baseline
util.cnc rf3_baseline
show cartoon, rf3_baseline

load ./rf3_6RP1_0.01.cif, rf3_gw_0_01
align rf3_gw_0_01, density_input
color rf3_gold, rf3_gw_0_01
util.cnc rf3_gw_0_01
show cartoon, rf3_gw_0_01

# Groups for easy show/hide toggling
group protenix, protenix_baseline protenix_gw_0_1
group rf3, rf3_baseline rf3_gw_0_01

# Hide chain-break dashes
set cartoon_gap_cutoff, 0

# Background / depth (override .pymolrc as needed)
bg_color white
set ray_shadow, 0
set ambient_occlusion_mode, 0

# Set all state
set state, 0, all

# ---------------------------------------------------------------
# VIEWING ANGLE
# ---------------------------------------------------------------


# ---------------------------------------------------------------
# RAY TRACE
# Toggle visibility per-render, e.g.:
hide everything
show cartoon, protenix_baseline or density_input
ray 2400, 1800
png protenix_6RP1_baseline.png, dpi=300

hide everything, protenix_baseline
show cartoon, protenix_gw_0_1
ray 2400, 1800
png protenix_6RP1_guided.png, dpi=300

hide everything, protenix
show cartoon, rf3_baseline 
ray 2400, 1800
png rf3_6RP1_baseline.png, dpi=300

hide everything, rf3_baseline
show cartoon, rf3_gw_0_01
ray 2400, 1800
png rf3_6RP1_guided0_01.png, dpi=300
# ---------------------------------------------------------------

