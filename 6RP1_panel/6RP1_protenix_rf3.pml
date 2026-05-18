# 2QCV: density input + protenix + RF3 predictions
# Render protenix and rf3 separately by toggling show/hide before each ray.

# Custom color for rf3 (matches MODEL_PALETTE in style.py)
set_color rf3_gold, [1.0, 0.843, 0.0]

# Density input (reference)
load ./2QCV_single_001_density_input.cif, density_input
color gray70, density_input
util.cnc density_input
show cartoon, density_input

# Protenix predictions
load ./protenix_2QCV_0.0.cif, protenix_baseline
align protenix_baseline, density_input
color good_pink, protenix_baseline
util.cnc protenix_baseline
show cartoon, protenix_baseline

load ./protenix_2QCV_0.1.cif, protenix_gw_0_1
align protenix_gw_0_1, density_input
color good_pink, protenix_gw_0_1
util.cnc protenix_gw_0_1
show cartoon, protenix_gw_0_1

# RF3 predictions
load ./rf3_2QCV_0.0.cif, rf3_baseline
align rf3_baseline, density_input
color rf3_gold, rf3_baseline
util.cnc rf3_baseline
show cartoon, rf3_baseline

load ./rf3_2QCV_0.005.cif, rf3_gw_0_005
align rf3_gw_0_005, density_input
color rf3_gold, rf3_gw_0_005
util.cnc rf3_gw_0_005
show cartoon, rf3_gw_0_005

load ./rf3_2QCV_0.01.cif, rf3_gw_0_01
align rf3_gw_0_01, density_input
color rf3_gold, rf3_gw_0_01
util.cnc rf3_gw_0_01
show cartoon, rf3_gw_0_01

load ./rf3_2QCV_0.02.cif, rf3_gw_0_02
align rf3_gw_0_02, density_input
color rf3_gold, rf3_gw_0_02
util.cnc rf3_gw_0_02
show cartoon, rf3_gw_0_02

# Groups for easy show/hide toggling
group protenix, protenix_baseline protenix_gw_0_1
group rf3, rf3_baseline rf3_gw_0_005 rf3_gw_0_01 rf3_gw_0_02

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
set_view (\
    -0.361618608,    0.931795418,   -0.031427968,\
     0.524363399,    0.231139973,    0.819521427,\
     0.770892859,    0.279874563,   -0.572183430,\
     0.000000000,    0.000000000,  -64.002029419,\
    32.452522278,   51.082225800,   14.793140411,\
    55.797588348,   72.206466675,   20.000000000 )

# ---------------------------------------------------------------
# RAY TRACE
# Toggle visibility per-render, e.g.:
hide everything
show cartoon, protenix_baseline or density_input
ray 2400, 1800
png protenix_2QCV_baseline.png, dpi=300

hide everything, protenix_baseline
show cartoon, protenix_gw_0_1
ray 2400, 1800
png protenix_2QCV_guided.png, dpi=300

hide everything, protenix
show cartoon, rf3_baseline 
ray 2400, 1800
png rf3_2QCV_baseline.png, dpi=300

hide everything, rf3_baseline
show cartoon, rf3_gw_0_02
ray 2400, 1800
png rf3_2QCV_guided0_02.png, dpi=300
# ---------------------------------------------------------------

