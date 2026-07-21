# 6NI6 0.5occA_0.5occB panel - residue 201-206 (chain A) focus region
# 4 MAIN models only - boltz2 X-ray, boltz2 MD, protenix, rf3
# Each model shown baseline (gw0.0) vs guided
#   boltz2 / protenix : gw0.0  gw0.1
#   rf3               : gw0.0  gw0.01
# resi 201-206 (chain A) is a small-loop alt-loc region (alt A/B at 201-202 and
# 204-206 per classify_altloc_regions.csv). Change RESI_FOCUS below to retarget.
# Each model object is an ens8 stack (8 states) - baseline and guided share one
# color - distinguish them by toggling show/hide before each ray

# ---------------------------------------------------------------
# Custom colors - mirror ~/.pymolrc so this script is self-contained
# ---------------------------------------------------------------
set_color good_darkblue, [75, 95, 170]
set_color good_teal,     [0.310, 0.725, 0.686]
set_color good_pink,     [255, 172, 183]
set_color good_navaho,   [255, 224, 172]

# ---------------------------------------------------------------
# Alt-loc cartoon fix (non-rebond default for PyMOL < 3.2)
# These crystallographic density inputs encode the blank alt-id as a quoted
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

# ---------------------------------------------------------------
# Density input (reference) - stays in the map frame, keeps alt A and alt B
# ---------------------------------------------------------------
load_altfix ./6NI6_density_input.cif, density_input
color gray70, density_input
util.cnc density_input
show cartoon, density_input

# ---------------------------------------------------------------
# Boltz2 X-ray density   (good_darkblue)
# ---------------------------------------------------------------
load ./6NI6_boltz2_xray_gw0p0.cif, boltz2_xray_baseline
align boltz2_xray_baseline, density_input
color good_darkblue, boltz2_xray_baseline
util.cnc boltz2_xray_baseline
show cartoon, boltz2_xray_baseline

load ./6NI6_boltz2_xray_gw0p1.cif, boltz2_xray_guided
align boltz2_xray_guided, density_input
color good_darkblue, boltz2_xray_guided
util.cnc boltz2_xray_guided
show cartoon, boltz2_xray_guided

# ---------------------------------------------------------------
# Boltz2 MD density   (good_teal)
# ---------------------------------------------------------------
load ./6NI6_boltz2_md_gw0p0.cif, boltz2_md_baseline
align boltz2_md_baseline, density_input
color good_teal, boltz2_md_baseline
util.cnc boltz2_md_baseline
show cartoon, boltz2_md_baseline

load ./6NI6_boltz2_md_gw0p1.cif, boltz2_md_guided
align boltz2_md_guided, density_input
color good_teal, boltz2_md_guided
util.cnc boltz2_md_guided
show cartoon, boltz2_md_guided

# ---------------------------------------------------------------
# Protenix   (good_pink)
# ---------------------------------------------------------------
load ./6NI6_protenix_gw0p0.cif, protenix_baseline
align protenix_baseline, density_input
color good_pink, protenix_baseline
util.cnc protenix_baseline
show cartoon, protenix_baseline

load ./6NI6_protenix_gw0p1.cif, protenix_guided
align protenix_guided, density_input
color good_pink, protenix_guided
util.cnc protenix_guided
show cartoon, protenix_guided

# ---------------------------------------------------------------
# RF3 fixed   (good_navaho)
# ---------------------------------------------------------------
load ./6NI6_rf3_fixed_gw0p0.cif, rf3_fixed_baseline
align rf3_fixed_baseline, density_input
color good_navaho, rf3_fixed_baseline
util.cnc rf3_fixed_baseline
show cartoon, rf3_fixed_baseline

load ./6NI6_rf3_fixed_gw0p01.cif, rf3_fixed_guided
align rf3_fixed_guided, density_input
color good_navaho, rf3_fixed_guided
util.cnc rf3_fixed_guided
show cartoon, rf3_fixed_guided

# ---------------------------------------------------------------
# Groups for easy show/hide toggling (flat - one level only)
# ---------------------------------------------------------------
group boltz2_xray, boltz2_xray_baseline boltz2_xray_guided
group boltz2_md, boltz2_md_baseline boltz2_md_guided
group protenix, protenix_baseline protenix_guided
group rf3_fixed, rf3_fixed_baseline rf3_fixed_guided

# ---------------------------------------------------------------
# Display settings
# ---------------------------------------------------------------
# Hide chain-break dashes
set cartoon_gap_cutoff, 0
# Render all model cartoons as loop - removes helix and sheet secondary structure
cartoon loop, not density_input
# Background and depth - override .pymolrc as needed
bg_color white
set ray_shadow, 0
set ambient_occlusion_mode, 0
# Cartoon outlines (BLOPIG method) - normal color plus an edge line during ray
# ray_trace_mode 1 = color + outline. Higher ray_trace_gain = thicker/darker edge.
# Outline color is inherited from .pymolrc (gray30) - antialias 2 also from there.
# Note: transparent objects receive NO outline in mode 1, so only opaque objects
# are edged (see the transparency block below for which objects that is).
set ray_trace_mode, 1
set ray_trace_gain, 1.5
# Show all 8 ensemble states at once
set state, 0, all

# Per-residue focus - resi 201-207 (chain A) side-chain sticks on every object
show sticks, resi 201-207

# Transparency - models see-through, density_input stays opaque
set cartoon_transparency, 0.3
set cartoon_transparency, 0.0, density_input
set stick_transparency, 0.3
set stick_transparency, 0.0, density_input

# ---------------------------------------------------------------
# VIEWING ANGLE - FILL IN
# Temporary auto-framing on the residues of interest - replace with a pinned view.
orient chain A and resi 201-206
zoom chain A and resi 201-206, 4
# ---- PASTE YOUR set_view HERE -----------------------------------
# Fine-tune the angle in the GUI, run  get_view , then paste the 18-number
# matrix below and uncomment it. A captured set_view supersedes the orient/zoom
# above, so it must stay AFTER them.

set_view (\
     0.899482489,   -0.293715417,    0.323512465,\
     0.355351508,    0.922548532,   -0.150431931,\
    -0.254271835,    0.250271916,    0.934189200,\
     0.000000000,    0.000000000,  -47.149951935,\
    44.778808594,   16.593805313,   43.532207489,\
    38.437721252,   55.862182617,   20.000000000 )

# ---------------------------------------------------------------
# RAY TRACE / EXPORT
# Renders one model condition at a time over the gray reference
# Density input stays visible throughout
# Comment out any family block you do not need - each ray is full-res
# ---------------------------------------------------------------
disable all
enable density_input

# Un-mask the group containers - a disabled parent group hides its members at
# render time even when the member object itself is enabled. disable all turns
# the groups off too, so enabling a leaf below would never render without this.
# The leaves stay individually disabled (from disable all) until toggled per block.
enable boltz2_xray
enable boltz2_md
enable protenix
enable rf3_fixed

# --- Boltz2 X-ray  (baseline gw0.0 -> guided gw0.1) ---
enable boltz2_xray_baseline
ray 2400, 1800
png 6NI6_boltz2_xray_baseline.png, dpi=300
disable boltz2_xray_baseline
enable boltz2_xray_guided
ray 2400, 1800
png 6NI6_boltz2_xray_guided.png, dpi=300
disable boltz2_xray_guided

# --- Boltz2 MD  (baseline gw0.0 -> guided gw0.1) ---
enable boltz2_md_baseline
ray 2400, 1800
png 6NI6_boltz2_md_baseline.png, dpi=300
disable boltz2_md_baseline
enable boltz2_md_guided
ray 2400, 1800
png 6NI6_boltz2_md_guided.png, dpi=300
disable boltz2_md_guided

# --- Protenix  (baseline gw0.0 -> guided gw0.1) ---
enable protenix_baseline
ray 2400, 1800
png 6NI6_protenix_baseline.png, dpi=300
disable protenix_baseline
enable protenix_guided
ray 2400, 1800
png 6NI6_protenix_guided.png, dpi=300
disable protenix_guided

# --- RF3 fixed  (baseline gw0.0 -> guided gw0.01) ---
enable rf3_fixed_baseline
ray 2400, 1800
png 6NI6_rf3_fixed_baseline.png, dpi=300
disable rf3_fixed_baseline
enable rf3_fixed_guided
ray 2400, 1800
png 6NI6_rf3_fixed_guided.png, dpi=300
disable rf3_fixed_guided

# Re-enable everything for interactive inspection
enable all
