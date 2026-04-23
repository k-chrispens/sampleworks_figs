# Load baseline predictions for 6B8X and 2A26, aligned to ground truth

# ============================================================
# Ground truth
# ============================================================
load ../../../processed/6B8X/6B8X_single_001_density_input.cif, 6B8X_gt
color good_gray, 6B8X_gt
util.cnc 6B8X_gt
show cartoon, 6B8X_gt

load ../../../processed/2A26/2A26_single_001_density_input.cif, 2A26_gt
# Fix altloc-induced chain breaks: resi 12(' ')-13(A/B) and resi 20(A/B)-21(' ')
bond 2A26_gt and resi 12 and name C and not alt B, 2A26_gt and resi 13 and name N and alt A
bond 2A26_gt and resi 12 and name C and not alt A, 2A26_gt and resi 13 and name N and alt B
bond 2A26_gt and resi 20 and name C and alt A, 2A26_gt and resi 21 and name N and not alt B
bond 2A26_gt and resi 20 and name C and alt B, 2A26_gt and resi 21 and name N and not alt A
color good_gray, 2A26_gt
util.cnc 2A26_gt
show cartoon, 2A26_gt

# ============================================================
# 6B8X predictions
# ============================================================
# --- 6B8X boltz2_XRAY (good_darkblue) ---
load ./boltz2_XRAY/6B8X_native_occ/seed_42/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2XRAY_s42
align 6B8X_boltz2XRAY_s42, 6B8X_gt
color good_darkblue, 6B8X_boltz2XRAY_s42
util.cnc 6B8X_boltz2XRAY_s42
show cartoon, 6B8X_boltz2XRAY_s42
load ./boltz2_XRAY/6B8X_native_occ/seed_43/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2XRAY_s43
align 6B8X_boltz2XRAY_s43, 6B8X_gt
color good_darkblue, 6B8X_boltz2XRAY_s43
util.cnc 6B8X_boltz2XRAY_s43
show cartoon, 6B8X_boltz2XRAY_s43
load ./boltz2_XRAY/6B8X_native_occ/seed_44/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2XRAY_s44
align 6B8X_boltz2XRAY_s44, 6B8X_gt
color good_darkblue, 6B8X_boltz2XRAY_s44
util.cnc 6B8X_boltz2XRAY_s44
show cartoon, 6B8X_boltz2XRAY_s44
load ./boltz2_XRAY/6B8X_native_occ/seed_45/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2XRAY_s45
align 6B8X_boltz2XRAY_s45, 6B8X_gt
color good_darkblue, 6B8X_boltz2XRAY_s45
util.cnc 6B8X_boltz2XRAY_s45
show cartoon, 6B8X_boltz2XRAY_s45
load ./boltz2_XRAY/6B8X_native_occ/seed_46/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2XRAY_s46
align 6B8X_boltz2XRAY_s46, 6B8X_gt
color good_darkblue, 6B8X_boltz2XRAY_s46
util.cnc 6B8X_boltz2XRAY_s46
show cartoon, 6B8X_boltz2XRAY_s46
load ./boltz2_XRAY/6B8X_native_occ/seed_47/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2XRAY_s47
align 6B8X_boltz2XRAY_s47, 6B8X_gt
color good_darkblue, 6B8X_boltz2XRAY_s47
util.cnc 6B8X_boltz2XRAY_s47
show cartoon, 6B8X_boltz2XRAY_s47
load ./boltz2_XRAY/6B8X_native_occ/seed_48/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2XRAY_s48
align 6B8X_boltz2XRAY_s48, 6B8X_gt
color good_darkblue, 6B8X_boltz2XRAY_s48
util.cnc 6B8X_boltz2XRAY_s48
show cartoon, 6B8X_boltz2XRAY_s48
load ./boltz2_XRAY/6B8X_native_occ/seed_49/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2XRAY_s49
align 6B8X_boltz2XRAY_s49, 6B8X_gt
color good_darkblue, 6B8X_boltz2XRAY_s49
util.cnc 6B8X_boltz2XRAY_s49
show cartoon, 6B8X_boltz2XRAY_s49

# --- 6B8X boltz2_MD (good_teal) ---
load ./boltz2_MD/6B8X_native_occ/seed_42/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2MD_s42
align 6B8X_boltz2MD_s42, 6B8X_gt
color good_teal, 6B8X_boltz2MD_s42
util.cnc 6B8X_boltz2MD_s42
show cartoon, 6B8X_boltz2MD_s42
load ./boltz2_MD/6B8X_native_occ/seed_43/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2MD_s43
align 6B8X_boltz2MD_s43, 6B8X_gt
color good_teal, 6B8X_boltz2MD_s43
util.cnc 6B8X_boltz2MD_s43
show cartoon, 6B8X_boltz2MD_s43
load ./boltz2_MD/6B8X_native_occ/seed_44/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2MD_s44
align 6B8X_boltz2MD_s44, 6B8X_gt
color good_teal, 6B8X_boltz2MD_s44
util.cnc 6B8X_boltz2MD_s44
show cartoon, 6B8X_boltz2MD_s44
load ./boltz2_MD/6B8X_native_occ/seed_45/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2MD_s45
align 6B8X_boltz2MD_s45, 6B8X_gt
color good_teal, 6B8X_boltz2MD_s45
util.cnc 6B8X_boltz2MD_s45
show cartoon, 6B8X_boltz2MD_s45
load ./boltz2_MD/6B8X_native_occ/seed_46/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2MD_s46
align 6B8X_boltz2MD_s46, 6B8X_gt
color good_teal, 6B8X_boltz2MD_s46
util.cnc 6B8X_boltz2MD_s46
show cartoon, 6B8X_boltz2MD_s46
load ./boltz2_MD/6B8X_native_occ/seed_47/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2MD_s47
align 6B8X_boltz2MD_s47, 6B8X_gt
color good_teal, 6B8X_boltz2MD_s47
util.cnc 6B8X_boltz2MD_s47
show cartoon, 6B8X_boltz2MD_s47
load ./boltz2_MD/6B8X_native_occ/seed_48/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2MD_s48
align 6B8X_boltz2MD_s48, 6B8X_gt
color good_teal, 6B8X_boltz2MD_s48
util.cnc 6B8X_boltz2MD_s48
show cartoon, 6B8X_boltz2MD_s48
load ./boltz2_MD/6B8X_native_occ/seed_49/boltz_results_6B8X_native_occ/predictions/6B8X_native_occ/6B8X_native_occ_model_0.cif, 6B8X_boltz2MD_s49
align 6B8X_boltz2MD_s49, 6B8X_gt
color good_teal, 6B8X_boltz2MD_s49
util.cnc 6B8X_boltz2MD_s49
show cartoon, 6B8X_boltz2MD_s49

# --- 6B8X rf3 (good_navaho) ---
load ./rf3/6B8X_native_occ/seed_42/6B8X_single_001_density_input/seed-42_sample-0/6B8X_single_001_density_input_seed-42_sample-0_model.cif, 6B8X_rf3_s42
align 6B8X_rf3_s42, 6B8X_gt
color good_navaho, 6B8X_rf3_s42
util.cnc 6B8X_rf3_s42
show cartoon, 6B8X_rf3_s42
load ./rf3/6B8X_native_occ/seed_43/6B8X_single_001_density_input/seed-43_sample-0/6B8X_single_001_density_input_seed-43_sample-0_model.cif, 6B8X_rf3_s43
align 6B8X_rf3_s43, 6B8X_gt
color good_navaho, 6B8X_rf3_s43
util.cnc 6B8X_rf3_s43
show cartoon, 6B8X_rf3_s43
load ./rf3/6B8X_native_occ/seed_44/6B8X_single_001_density_input/seed-44_sample-0/6B8X_single_001_density_input_seed-44_sample-0_model.cif, 6B8X_rf3_s44
align 6B8X_rf3_s44, 6B8X_gt
color good_navaho, 6B8X_rf3_s44
util.cnc 6B8X_rf3_s44
show cartoon, 6B8X_rf3_s44
load ./rf3/6B8X_native_occ/seed_45/6B8X_single_001_density_input/seed-45_sample-0/6B8X_single_001_density_input_seed-45_sample-0_model.cif, 6B8X_rf3_s45
align 6B8X_rf3_s45, 6B8X_gt
color good_navaho, 6B8X_rf3_s45
util.cnc 6B8X_rf3_s45
show cartoon, 6B8X_rf3_s45
load ./rf3/6B8X_native_occ/seed_46/6B8X_single_001_density_input/seed-46_sample-0/6B8X_single_001_density_input_seed-46_sample-0_model.cif, 6B8X_rf3_s46
align 6B8X_rf3_s46, 6B8X_gt
color good_navaho, 6B8X_rf3_s46
util.cnc 6B8X_rf3_s46
show cartoon, 6B8X_rf3_s46
load ./rf3/6B8X_native_occ/seed_47/6B8X_single_001_density_input/seed-47_sample-0/6B8X_single_001_density_input_seed-47_sample-0_model.cif, 6B8X_rf3_s47
align 6B8X_rf3_s47, 6B8X_gt
color good_navaho, 6B8X_rf3_s47
util.cnc 6B8X_rf3_s47
show cartoon, 6B8X_rf3_s47
load ./rf3/6B8X_native_occ/seed_48/6B8X_single_001_density_input/seed-48_sample-0/6B8X_single_001_density_input_seed-48_sample-0_model.cif, 6B8X_rf3_s48
align 6B8X_rf3_s48, 6B8X_gt
color good_navaho, 6B8X_rf3_s48
util.cnc 6B8X_rf3_s48
show cartoon, 6B8X_rf3_s48
load ./rf3/6B8X_native_occ/seed_49/6B8X_single_001_density_input/seed-49_sample-0/6B8X_single_001_density_input_seed-49_sample-0_model.cif, 6B8X_rf3_s49
align 6B8X_rf3_s49, 6B8X_gt
color good_navaho, 6B8X_rf3_s49
util.cnc 6B8X_rf3_s49
show cartoon, 6B8X_rf3_s49

# --- 6B8X protenix (good_pink) ---
load ./protenix/6B8X_native_occ/6B8X_native_occ/seed_42/predictions/6B8X_native_occ_sample_0.cif, 6B8X_protenix_s42
align 6B8X_protenix_s42, 6B8X_gt
color good_pink, 6B8X_protenix_s42
util.cnc 6B8X_protenix_s42
show cartoon, 6B8X_protenix_s42
load ./protenix/6B8X_native_occ/6B8X_native_occ/seed_43/predictions/6B8X_native_occ_sample_0.cif, 6B8X_protenix_s43
align 6B8X_protenix_s43, 6B8X_gt
color good_pink, 6B8X_protenix_s43
util.cnc 6B8X_protenix_s43
show cartoon, 6B8X_protenix_s43
load ./protenix/6B8X_native_occ/6B8X_native_occ/seed_44/predictions/6B8X_native_occ_sample_0.cif, 6B8X_protenix_s44
align 6B8X_protenix_s44, 6B8X_gt
color good_pink, 6B8X_protenix_s44
util.cnc 6B8X_protenix_s44
show cartoon, 6B8X_protenix_s44
load ./protenix/6B8X_native_occ/6B8X_native_occ/seed_45/predictions/6B8X_native_occ_sample_0.cif, 6B8X_protenix_s45
align 6B8X_protenix_s45, 6B8X_gt
color good_pink, 6B8X_protenix_s45
util.cnc 6B8X_protenix_s45
show cartoon, 6B8X_protenix_s45
load ./protenix/6B8X_native_occ/6B8X_native_occ/seed_46/predictions/6B8X_native_occ_sample_0.cif, 6B8X_protenix_s46
align 6B8X_protenix_s46, 6B8X_gt
color good_pink, 6B8X_protenix_s46
util.cnc 6B8X_protenix_s46
show cartoon, 6B8X_protenix_s46
load ./protenix/6B8X_native_occ/6B8X_native_occ/seed_47/predictions/6B8X_native_occ_sample_0.cif, 6B8X_protenix_s47
align 6B8X_protenix_s47, 6B8X_gt
color good_pink, 6B8X_protenix_s47
util.cnc 6B8X_protenix_s47
show cartoon, 6B8X_protenix_s47
load ./protenix/6B8X_native_occ/6B8X_native_occ/seed_48/predictions/6B8X_native_occ_sample_0.cif, 6B8X_protenix_s48
align 6B8X_protenix_s48, 6B8X_gt
color good_pink, 6B8X_protenix_s48
util.cnc 6B8X_protenix_s48
show cartoon, 6B8X_protenix_s48
load ./protenix/6B8X_native_occ/6B8X_native_occ/seed_49/predictions/6B8X_native_occ_sample_0.cif, 6B8X_protenix_s49
align 6B8X_protenix_s49, 6B8X_gt
color good_pink, 6B8X_protenix_s49
util.cnc 6B8X_protenix_s49
show cartoon, 6B8X_protenix_s49

# ============================================================
# 2A26 predictions
# ============================================================
# --- 2A26 boltz2_XRAY (good_darkblue) ---
load ./boltz2_XRAY/2A26_native_occ/seed_42/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2XRAY_s42
align 2A26_boltz2XRAY_s42, 2A26_gt
color good_darkblue, 2A26_boltz2XRAY_s42
util.cnc 2A26_boltz2XRAY_s42
show cartoon, 2A26_boltz2XRAY_s42
load ./boltz2_XRAY/2A26_native_occ/seed_43/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2XRAY_s43
align 2A26_boltz2XRAY_s43, 2A26_gt
color good_darkblue, 2A26_boltz2XRAY_s43
util.cnc 2A26_boltz2XRAY_s43
show cartoon, 2A26_boltz2XRAY_s43
load ./boltz2_XRAY/2A26_native_occ/seed_44/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2XRAY_s44
align 2A26_boltz2XRAY_s44, 2A26_gt
color good_darkblue, 2A26_boltz2XRAY_s44
util.cnc 2A26_boltz2XRAY_s44
show cartoon, 2A26_boltz2XRAY_s44
load ./boltz2_XRAY/2A26_native_occ/seed_45/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2XRAY_s45
align 2A26_boltz2XRAY_s45, 2A26_gt
color good_darkblue, 2A26_boltz2XRAY_s45
util.cnc 2A26_boltz2XRAY_s45
show cartoon, 2A26_boltz2XRAY_s45
load ./boltz2_XRAY/2A26_native_occ/seed_46/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2XRAY_s46
align 2A26_boltz2XRAY_s46, 2A26_gt
color good_darkblue, 2A26_boltz2XRAY_s46
util.cnc 2A26_boltz2XRAY_s46
show cartoon, 2A26_boltz2XRAY_s46
load ./boltz2_XRAY/2A26_native_occ/seed_47/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2XRAY_s47
align 2A26_boltz2XRAY_s47, 2A26_gt
color good_darkblue, 2A26_boltz2XRAY_s47
util.cnc 2A26_boltz2XRAY_s47
show cartoon, 2A26_boltz2XRAY_s47
load ./boltz2_XRAY/2A26_native_occ/seed_48/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2XRAY_s48
align 2A26_boltz2XRAY_s48, 2A26_gt
color good_darkblue, 2A26_boltz2XRAY_s48
util.cnc 2A26_boltz2XRAY_s48
show cartoon, 2A26_boltz2XRAY_s48
load ./boltz2_XRAY/2A26_native_occ/seed_49/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2XRAY_s49
align 2A26_boltz2XRAY_s49, 2A26_gt
color good_darkblue, 2A26_boltz2XRAY_s49
util.cnc 2A26_boltz2XRAY_s49
show cartoon, 2A26_boltz2XRAY_s49

# --- 2A26 boltz2_MD (good_teal) ---
load ./boltz2_MD/2A26_native_occ/seed_42/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2MD_s42
align 2A26_boltz2MD_s42, 2A26_gt
color good_teal, 2A26_boltz2MD_s42
util.cnc 2A26_boltz2MD_s42
show cartoon, 2A26_boltz2MD_s42
load ./boltz2_MD/2A26_native_occ/seed_43/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2MD_s43
align 2A26_boltz2MD_s43, 2A26_gt
color good_teal, 2A26_boltz2MD_s43
util.cnc 2A26_boltz2MD_s43
show cartoon, 2A26_boltz2MD_s43
load ./boltz2_MD/2A26_native_occ/seed_44/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2MD_s44
align 2A26_boltz2MD_s44, 2A26_gt
color good_teal, 2A26_boltz2MD_s44
util.cnc 2A26_boltz2MD_s44
show cartoon, 2A26_boltz2MD_s44
load ./boltz2_MD/2A26_native_occ/seed_45/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2MD_s45
align 2A26_boltz2MD_s45, 2A26_gt
color good_teal, 2A26_boltz2MD_s45
util.cnc 2A26_boltz2MD_s45
show cartoon, 2A26_boltz2MD_s45
load ./boltz2_MD/2A26_native_occ/seed_46/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2MD_s46
align 2A26_boltz2MD_s46, 2A26_gt
color good_teal, 2A26_boltz2MD_s46
util.cnc 2A26_boltz2MD_s46
show cartoon, 2A26_boltz2MD_s46
load ./boltz2_MD/2A26_native_occ/seed_47/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2MD_s47
align 2A26_boltz2MD_s47, 2A26_gt
color good_teal, 2A26_boltz2MD_s47
util.cnc 2A26_boltz2MD_s47
show cartoon, 2A26_boltz2MD_s47
load ./boltz2_MD/2A26_native_occ/seed_48/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2MD_s48
align 2A26_boltz2MD_s48, 2A26_gt
color good_teal, 2A26_boltz2MD_s48
util.cnc 2A26_boltz2MD_s48
show cartoon, 2A26_boltz2MD_s48
load ./boltz2_MD/2A26_native_occ/seed_49/boltz_results_2A26_native_occ/predictions/2A26_native_occ/2A26_native_occ_model_0.cif, 2A26_boltz2MD_s49
align 2A26_boltz2MD_s49, 2A26_gt
color good_teal, 2A26_boltz2MD_s49
util.cnc 2A26_boltz2MD_s49
show cartoon, 2A26_boltz2MD_s49

# --- 2A26 rf3 (good_navaho) ---
load ./rf3/2A26_native_occ/seed_42/2A26_single_001_density_input/seed-42_sample-0/2A26_single_001_density_input_seed-42_sample-0_model.cif, 2A26_rf3_s42
align 2A26_rf3_s42, 2A26_gt
color good_navaho, 2A26_rf3_s42
util.cnc 2A26_rf3_s42
show cartoon, 2A26_rf3_s42
load ./rf3/2A26_native_occ/seed_43/2A26_single_001_density_input/seed-43_sample-0/2A26_single_001_density_input_seed-43_sample-0_model.cif, 2A26_rf3_s43
align 2A26_rf3_s43, 2A26_gt
color good_navaho, 2A26_rf3_s43
util.cnc 2A26_rf3_s43
show cartoon, 2A26_rf3_s43
load ./rf3/2A26_native_occ/seed_44/2A26_single_001_density_input/seed-44_sample-0/2A26_single_001_density_input_seed-44_sample-0_model.cif, 2A26_rf3_s44
align 2A26_rf3_s44, 2A26_gt
color good_navaho, 2A26_rf3_s44
util.cnc 2A26_rf3_s44
show cartoon, 2A26_rf3_s44
load ./rf3/2A26_native_occ/seed_45/2A26_single_001_density_input/seed-45_sample-0/2A26_single_001_density_input_seed-45_sample-0_model.cif, 2A26_rf3_s45
align 2A26_rf3_s45, 2A26_gt
color good_navaho, 2A26_rf3_s45
util.cnc 2A26_rf3_s45
show cartoon, 2A26_rf3_s45
load ./rf3/2A26_native_occ/seed_46/2A26_single_001_density_input/seed-46_sample-0/2A26_single_001_density_input_seed-46_sample-0_model.cif, 2A26_rf3_s46
align 2A26_rf3_s46, 2A26_gt
color good_navaho, 2A26_rf3_s46
util.cnc 2A26_rf3_s46
show cartoon, 2A26_rf3_s46
load ./rf3/2A26_native_occ/seed_47/2A26_single_001_density_input/seed-47_sample-0/2A26_single_001_density_input_seed-47_sample-0_model.cif, 2A26_rf3_s47
align 2A26_rf3_s47, 2A26_gt
color good_navaho, 2A26_rf3_s47
util.cnc 2A26_rf3_s47
show cartoon, 2A26_rf3_s47
load ./rf3/2A26_native_occ/seed_48/2A26_single_001_density_input/seed-48_sample-0/2A26_single_001_density_input_seed-48_sample-0_model.cif, 2A26_rf3_s48
align 2A26_rf3_s48, 2A26_gt
color good_navaho, 2A26_rf3_s48
util.cnc 2A26_rf3_s48
show cartoon, 2A26_rf3_s48
load ./rf3/2A26_native_occ/seed_49/2A26_single_001_density_input/seed-49_sample-0/2A26_single_001_density_input_seed-49_sample-0_model.cif, 2A26_rf3_s49
align 2A26_rf3_s49, 2A26_gt
color good_navaho, 2A26_rf3_s49
util.cnc 2A26_rf3_s49
show cartoon, 2A26_rf3_s49

# --- 2A26 protenix (good_pink) ---
load ./protenix/2A26_native_occ/2A26_native_occ/seed_42/predictions/2A26_native_occ_sample_0.cif, 2A26_protenix_s42
align 2A26_protenix_s42, 2A26_gt
color good_pink, 2A26_protenix_s42
util.cnc 2A26_protenix_s42
show cartoon, 2A26_protenix_s42
load ./protenix/2A26_native_occ/2A26_native_occ/seed_43/predictions/2A26_native_occ_sample_0.cif, 2A26_protenix_s43
align 2A26_protenix_s43, 2A26_gt
color good_pink, 2A26_protenix_s43
util.cnc 2A26_protenix_s43
show cartoon, 2A26_protenix_s43
load ./protenix/2A26_native_occ/2A26_native_occ/seed_44/predictions/2A26_native_occ_sample_0.cif, 2A26_protenix_s44
align 2A26_protenix_s44, 2A26_gt
color good_pink, 2A26_protenix_s44
util.cnc 2A26_protenix_s44
show cartoon, 2A26_protenix_s44
load ./protenix/2A26_native_occ/2A26_native_occ/seed_45/predictions/2A26_native_occ_sample_0.cif, 2A26_protenix_s45
align 2A26_protenix_s45, 2A26_gt
color good_pink, 2A26_protenix_s45
util.cnc 2A26_protenix_s45
show cartoon, 2A26_protenix_s45
load ./protenix/2A26_native_occ/2A26_native_occ/seed_46/predictions/2A26_native_occ_sample_0.cif, 2A26_protenix_s46
align 2A26_protenix_s46, 2A26_gt
color good_pink, 2A26_protenix_s46
util.cnc 2A26_protenix_s46
show cartoon, 2A26_protenix_s46
load ./protenix/2A26_native_occ/2A26_native_occ/seed_47/predictions/2A26_native_occ_sample_0.cif, 2A26_protenix_s47
align 2A26_protenix_s47, 2A26_gt
color good_pink, 2A26_protenix_s47
util.cnc 2A26_protenix_s47
show cartoon, 2A26_protenix_s47
load ./protenix/2A26_native_occ/2A26_native_occ/seed_48/predictions/2A26_native_occ_sample_0.cif, 2A26_protenix_s48
align 2A26_protenix_s48, 2A26_gt
color good_pink, 2A26_protenix_s48
util.cnc 2A26_protenix_s48
show cartoon, 2A26_protenix_s48
load ./protenix/2A26_native_occ/2A26_native_occ/seed_49/predictions/2A26_native_occ_sample_0.cif, 2A26_protenix_s49
align 2A26_protenix_s49, 2A26_gt
color good_pink, 2A26_protenix_s49
util.cnc 2A26_protenix_s49
show cartoon, 2A26_protenix_s49

# ============================================================
# Groups — leaf groups (seeds into model groups)
# ============================================================
# 6B8X
group 6B8X_boltz2_XRAY, 6B8X_boltz2XRAY_s42 6B8X_boltz2XRAY_s43 6B8X_boltz2XRAY_s44 6B8X_boltz2XRAY_s45 6B8X_boltz2XRAY_s46 6B8X_boltz2XRAY_s47 6B8X_boltz2XRAY_s48 6B8X_boltz2XRAY_s49
group 6B8X_boltz2_MD, 6B8X_boltz2MD_s42 6B8X_boltz2MD_s43 6B8X_boltz2MD_s44 6B8X_boltz2MD_s45 6B8X_boltz2MD_s46 6B8X_boltz2MD_s47 6B8X_boltz2MD_s48 6B8X_boltz2MD_s49
group 6B8X_rf3, 6B8X_rf3_s42 6B8X_rf3_s43 6B8X_rf3_s44 6B8X_rf3_s45 6B8X_rf3_s46 6B8X_rf3_s47 6B8X_rf3_s48 6B8X_rf3_s49
group 6B8X_protenix, 6B8X_protenix_s42 6B8X_protenix_s43 6B8X_protenix_s44 6B8X_protenix_s45 6B8X_protenix_s46 6B8X_protenix_s47 6B8X_protenix_s48 6B8X_protenix_s49
# 2A26
group 2A26_boltz2_XRAY, 2A26_boltz2XRAY_s42 2A26_boltz2XRAY_s43 2A26_boltz2XRAY_s44 2A26_boltz2XRAY_s45 2A26_boltz2XRAY_s46 2A26_boltz2XRAY_s47 2A26_boltz2XRAY_s48 2A26_boltz2XRAY_s49
group 2A26_boltz2_MD, 2A26_boltz2MD_s42 2A26_boltz2MD_s43 2A26_boltz2MD_s44 2A26_boltz2MD_s45 2A26_boltz2MD_s46 2A26_boltz2MD_s47 2A26_boltz2MD_s48 2A26_boltz2MD_s49
group 2A26_rf3, 2A26_rf3_s42 2A26_rf3_s43 2A26_rf3_s44 2A26_rf3_s45 2A26_rf3_s46 2A26_rf3_s47 2A26_rf3_s48 2A26_rf3_s49
group 2A26_protenix, 2A26_protenix_s42 2A26_protenix_s43 2A26_protenix_s44 2A26_protenix_s45 2A26_protenix_s46 2A26_protenix_s47 2A26_protenix_s48 2A26_protenix_s49

# ============================================================
# Display settings
# ============================================================
set cartoon_transparency, 0.3
set cartoon_transparency, 0.0, *_gt
set ray_shadow, 0
# Black outlines on opaque GT via ray_trace_mode 1 (BLOPIG method)
# Transparent baselines will NOT get outlines in this mode
set ray_trace_mode, 1
set ray_trace_gain, 0.1
set antialias, 2
set cartoon_gap_cutoff, 3

# ============================================================
# 6B8X viewing
# ============================================================ 
set_view (\
     0.824064076,    0.132588267,   -0.550761402,\
    -0.438142091,    0.765453100,   -0.471287370,\
     0.359095037,    0.629683197,    0.688875377,\
    -0.000006092,   -0.000064261,  -42.700595856,\
    32.381988525,   44.253520966,   35.678207397,\
    26.406190872,   59.003044128,   20.000000000 )

ray 1680, 1262, async=1

# ============================================================
# 2A26 viewing
# ============================================================
show sticks, 2A26* and resi 14-20
set stick_transparency, 0.3, not *_gt

