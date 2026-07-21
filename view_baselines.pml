# Load baseline predictions for 6B8X and 2A26, aligned to ground truth

# ---------------------------------------------------------------
# Alt-loc cartoon fix (non-rebond default for PyMOL < 3.2)
# These crystallographic ground-truth inputs encode the blank alt-id as a
# quoted space (' ') instead of the mmCIF null (.). PyMOL only forms a
# blank<->A/B bond when the blank alt is the empty string, so every blank<->A/B
# backbone peptide bond is dropped at load and the cartoon shatters into
# fragments at alt-loc residues. The cartoon spline is cached from load-time
# connectivity, so a post-load bond cannot repair it - the file is fixed
# BEFORE load. load_altfix rewrites ONLY the label_alt_id column (quoted-space
# to dot) in memory, then loads, so blank<->A/B bonds form at load.
# PyMOL 3.2+ alternative (drop this block) is plain load, then
#   alter obj, alt=alt.strip()
#   rebond obj
# ---------------------------------------------------------------
python
def _altfix_text(txt):
    L = txt.split("\n")
    cols = [l.strip() for l in L if l.strip().startswith("_atom_site.")]
    if "_atom_site.label_alt_id" not in cols:
        return txt
    ai = cols.index("_atom_site.label_alt_id")
    def fix(line):
        i, n, k = 0, len(line), 0
        while i < n:
            while i < n and line[i] in " \t":
                i += 1
            if i >= n:
                break
            s = i
            if line[i] in "'\"":
                q = line[i]; i += 1
                while i < n and not (line[i] == q and (i + 1 >= n or line[i + 1] in " \t")):
                    i += 1
                i += 1; v = line[s + 1:i - 1]
            else:
                while i < n and line[i] not in " \t":
                    i += 1
                v = line[s:i]
            if k == ai:
                return (line[:s] + "." + line[i:]) if v.strip() == "" else line
            k += 1
        return line
    return "\n".join(fix(l) if (l[:4] == "ATOM" or l[:6] == "HETATM") else l for l in L)
def load_altfix(filename, oname):
    with open(filename) as _fh:
        cmd.load_raw(_altfix_text(_fh.read()), "cif", oname.strip())
cmd.extend("load_altfix", load_altfix)
python end

# ============================================================
# Ground truth
# ============================================================
load_altfix ../../../processed/6B8X/6B8X_single_001_density_input.cif, 6B8X_gt
color good_gray, 6B8X_gt
util.cnc 6B8X_gt
show cartoon, 6B8X_gt

load_altfix ../../../processed/2A26/2A26_single_001_density_input.cif, 2A26_gt
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
set ray_trace_gain, 1.5
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

