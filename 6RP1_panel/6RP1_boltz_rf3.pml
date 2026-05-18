# 6RP1: density input + boltzXRD + RF3 predictions
# Render boltzXRD and rf3 separately by toggling show/hide before each ray.

# Density input (reference)
load ./6RP1_single_001_density_input.cif, density_input
color gray70, density_input
util.cnc density_input
show cartoon, density_input

# boltzXRD predictions
load ./boltzXRD_6RP1_0.0.cif, boltzXRD_baseline
align boltzXRD_baseline, density_input
color good_pink, boltzXRD_baseline
util.cnc boltzXRD_baseline
show cartoon, boltzXRD_baseline

load ./boltzXRD_6RP1_0.1.cif, boltzXRD_gw_0_1
align boltzXRD_gw_0_1, density_input
color good_pink, boltzXRD_gw_0_1
util.cnc boltzXRD_gw_0_1
show cartoon, boltzXRD_gw_0_1

# RF3 predictions
load ./rf3_6RP1_0.0.cif, rf3_baseline
align rf3_baseline, density_input
color good_navaho, rf3_baseline
util.cnc rf3_baseline
show cartoon, rf3_baseline

load ./rf3_6RP1_0.01.cif, rf3_gw_0_01
align rf3_gw_0_01, density_input
color good_navaho, rf3_gw_0_01
util.cnc rf3_gw_0_01
show cartoon, rf3_gw_0_01

# Groups for easy show/hide toggling
group boltzXRD, boltzXRD_baseline boltzXRD_gw_0_1
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
# hide everything
# show cartoon, boltzXRD_baseline or density_input
# ray 2400, 1800
# png boltzXRD_6RP1_baseline.png, dpi=300
#
# hide everything, boltzXRD_baseline
# show cartoon, boltzXRD_gw_0_1
# ray 2400, 1800
# png boltzXRD_6RP1_guided.png, dpi=300
#
# hide everything, boltzXRD
# show cartoon, rf3_baseline
# ray 2400, 1800
# png rf3_6RP1_baseline.png, dpi=300
#
# hide everything, rf3_baseline
# show cartoon, rf3_gw_0_01
# ray 2400, 1800
# png rf3_6RP1_guided0_01.png, dpi=300
# ---------------------------------------------------------------
