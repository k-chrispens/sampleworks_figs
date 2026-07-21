# 2QCV: density input + boltz2 (X-ray, MD) + protenix + RF3 predictions
# Render each model separately by toggling show/hide before each ray.

# ---------------------------------------------------------------
# Alt-loc cartoon fix (non-rebond default for PyMOL < 3.2)
# This crystallographic density input encodes the blank alt-id as a quoted
# space (' ') instead of the mmCIF null (.). PyMOL only forms a blank<->A/B
# bond when the blank alt is the empty string, so every blank<->A/B backbone
# peptide bond is dropped at load and the cartoon shatters into fragments at
# alt-loc residues. The cartoon spline is cached from load-time connectivity,
# so a post-load bond cannot repair it - the file is fixed BEFORE load.
# load_altfix rewrites ONLY the label_alt_id column (quoted-space to dot) in
# memory, then loads, so blank<->A/B bonds form at load and the cartoon holds.
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

# Density input (reference)
load_altfix ./2QCV_single_001_density_input.cif, density_input
color gray70, density_input
util.cnc density_input
show cartoon, density_input

# Boltz2 X-ray predictions
load ./boltz2xray_2QCV_0.0.cif, boltz2_xray_baseline
align boltz2_xray_baseline, density_input
color good_darkblue, boltz2_xray_baseline
util.cnc boltz2_xray_baseline
show cartoon, boltz2_xray_baseline

load ./boltz2xray_2QCV_0.1.cif, boltz2_xray_gw_0_1
align boltz2_xray_gw_0_1, density_input
color good_darkblue, boltz2_xray_gw_0_1
util.cnc boltz2_xray_gw_0_1
show cartoon, boltz2_xray_gw_0_1

# Boltz2 MD predictions
load ./boltz2md_2QCV_0.0.cif, boltz2_md_baseline
align boltz2_md_baseline, density_input
color good_teal, boltz2_md_baseline
util.cnc boltz2_md_baseline
show cartoon, boltz2_md_baseline

load ./boltz2md_2QCV_0.1.cif, boltz2_md_gw_0_1
align boltz2_md_gw_0_1, density_input
color good_teal, boltz2_md_gw_0_1
util.cnc boltz2_md_gw_0_1
show cartoon, boltz2_md_gw_0_1

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
color good_navaho, rf3_baseline
util.cnc rf3_baseline
show cartoon, rf3_baseline

load ./rf3_2QCV_0.005.cif, rf3_gw_0_005
align rf3_gw_0_005, density_input
color good_navaho, rf3_gw_0_005
util.cnc rf3_gw_0_005
show cartoon, rf3_gw_0_005

load ./rf3_2QCV_0.01.cif, rf3_gw_0_01
align rf3_gw_0_01, density_input
color good_navaho, rf3_gw_0_01
util.cnc rf3_gw_0_01
show cartoon, rf3_gw_0_01

load ./rf3_2QCV_0.02.cif, rf3_gw_0_02
align rf3_gw_0_02, density_input
color good_navaho, rf3_gw_0_02
util.cnc rf3_gw_0_02
show cartoon, rf3_gw_0_02

# Groups for easy show/hide toggling
group boltz2_xray, boltz2_xray_baseline boltz2_xray_gw_0_1
group boltz2_md, boltz2_md_baseline boltz2_md_gw_0_1
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
disable all
enable density_input

# --- Boltz2 X-ray  (baseline gw0.0 -> guided gw0.1) ---
enable boltz2_xray
enable boltz2_xray_baseline
ray 2400, 1800
png boltz2_xray_2QCV_baseline.png, dpi=300

disable boltz2_xray_baseline
enable boltz2_xray_gw_0_1
ray 2400, 1800
png boltz2_xray_2QCV_guided.png, dpi=300

# --- Boltz2 MD  (baseline gw0.0 -> guided gw0.1) ---
disable boltz2_xray
enable boltz2_md
enable boltz2_md_baseline
ray 2400, 1800
png boltz2_md_2QCV_baseline.png, dpi=300

disable boltz2_md_baseline
enable boltz2_md_gw_0_1
ray 2400, 1800
png boltz2_md_2QCV_guided.png, dpi=300

# --- Protenix  (baseline gw0.0 -> guided gw0.1) ---
disable boltz2_md
enable protenix
enable protenix_baseline
ray 2400, 1800
png protenix_2QCV_baseline.png, dpi=300

disable protenix_baseline
enable protenix_gw_0_1
ray 2400, 1800
png protenix_2QCV_guided.png, dpi=300

disable protenix
enable rf3
enable rf3_baseline
ray 2400, 1800
png rf3_2QCV_baseline.png, dpi=300

disable rf3_baseline
enable rf3_gw_0_02
ray 2400, 1800
png rf3_2QCV_guided0_02.png, dpi=300

disable rf3_gw_0_02
enable rf3_gw_0_01
ray 2400, 1800
png rf3_2QCV_guided0_01.png, dpi=300
# ---------------------------------------------------------------
