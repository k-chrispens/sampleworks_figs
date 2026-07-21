# 5SOP 0.5occA_0.5occB panel - density input + boltz2 + protenix + rf3 ensembles
# Variant refined-patched.cif - each model object is an ens8 stack (8 states)
# Within a family the baseline (gw0.0) and guided share one color
# Distinguish baseline vs guided by toggling show/hide before each ray

# ---------------------------------------------------------------
# Custom colors - mirror ~/.pymolrc so this script is self-contained
# ---------------------------------------------------------------
set_color good_darkblue, [75, 95, 170]
set_color good_teal,     [0.310, 0.725, 0.686]
set_color good_pink,     [255, 172, 183]
set_color good_purple,   [213, 154, 181]
set_color good_melon,    [255, 198, 178]
set_color good_navaho,   [255, 224, 172]
# good_goldenrod is panel-specific (not in .pymolrc) - soft goldenrod so rf3-chiral-off reads distinct from rf3
set_color good_goldenrod, [227, 188, 88]

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
# Density input (reference) - stays in the map frame
# ---------------------------------------------------------------
load_altfix ./5SOP_density_input.cif, density_input
color gray70, density_input
remove density_input and alt C
util.cnc density_input
show cartoon, density_input

# ---------------------------------------------------------------
# Boltz2 X-ray density   (good_darkblue)
# ---------------------------------------------------------------
load ./5SOP_boltz2_xray_gw0p0.cif, boltz2_xray_baseline
align boltz2_xray_baseline, density_input
color good_darkblue, boltz2_xray_baseline
util.cnc boltz2_xray_baseline
show cartoon, boltz2_xray_baseline

load ./5SOP_boltz2_xray_gw0p1.cif, boltz2_xray_gw_0_1
align boltz2_xray_gw_0_1, density_input
color good_darkblue, boltz2_xray_gw_0_1
util.cnc boltz2_xray_gw_0_1
show cartoon, boltz2_xray_gw_0_1

# ---------------------------------------------------------------
# Boltz2 MD density   (good_teal)
# ---------------------------------------------------------------
load ./5SOP_boltz2_md_gw0p0.cif, boltz2_md_baseline
align boltz2_md_baseline, density_input
color good_teal, boltz2_md_baseline
util.cnc boltz2_md_baseline
show cartoon, boltz2_md_baseline

load ./5SOP_boltz2_md_gw0p1.cif, boltz2_md_gw_0_1
align boltz2_md_gw_0_1, density_input
color good_teal, boltz2_md_gw_0_1
util.cnc boltz2_md_gw_0_1
show cartoon, boltz2_md_gw_0_1

# ---------------------------------------------------------------
# Protenix   (good_pink)
# ---------------------------------------------------------------
load ./5SOP_protenix_gw0p0.cif, protenix_baseline
align protenix_baseline, density_input
color good_pink, protenix_baseline
util.cnc protenix_baseline
show cartoon, protenix_baseline

load ./5SOP_protenix_gw0p1.cif, protenix_gw_0_1
align protenix_gw_0_1, density_input
color good_pink, protenix_gw_0_1
util.cnc protenix_gw_0_1
show cartoon, protenix_gw_0_1

# ---------------------------------------------------------------
# Protenix mini   (good_purple)
# ---------------------------------------------------------------
load ./5SOP_protenix_mini_gw0p0.cif, protenix_mini_baseline
align protenix_mini_baseline, density_input
color good_purple, protenix_mini_baseline
util.cnc protenix_mini_baseline
show cartoon, protenix_mini_baseline

load ./5SOP_protenix_mini_gw0p1.cif, protenix_mini_gw_0_1
align protenix_mini_gw_0_1, density_input
color good_purple, protenix_mini_gw_0_1
util.cnc protenix_mini_gw_0_1
show cartoon, protenix_mini_gw_0_1

# ---------------------------------------------------------------
# Protenix tiny   (good_melon)
# ---------------------------------------------------------------
load ./5SOP_protenix_tiny_gw0p0.cif, protenix_tiny_baseline
align protenix_tiny_baseline, density_input
color good_melon, protenix_tiny_baseline
util.cnc protenix_tiny_baseline
show cartoon, protenix_tiny_baseline

load ./5SOP_protenix_tiny_gw0p1.cif, protenix_tiny_gw_0_1
align protenix_tiny_gw_0_1, density_input
color good_melon, protenix_tiny_gw_0_1
util.cnc protenix_tiny_gw_0_1
show cartoon, protenix_tiny_gw_0_1

# ---------------------------------------------------------------
# RF3 fixed   (good_navaho)
# ---------------------------------------------------------------
load ./5SOP_rf3_fixed_gw0p0.cif, rf3_fixed_baseline
align rf3_fixed_baseline, density_input
color good_navaho, rf3_fixed_baseline
util.cnc rf3_fixed_baseline
show cartoon, rf3_fixed_baseline

load ./5SOP_rf3_fixed_gw0p01.cif, rf3_fixed_gw_0_01
align rf3_fixed_gw_0_01, density_input
color good_navaho, rf3_fixed_gw_0_01
util.cnc rf3_fixed_gw_0_01
show cartoon, rf3_fixed_gw_0_01

# ---------------------------------------------------------------
# RF3 chiral-off   (good_goldenrod)
# ---------------------------------------------------------------
load ./5SOP_rf3_chiral_off_gw0p0.cif, rf3_chiral_off_baseline
align rf3_chiral_off_baseline, density_input
color good_goldenrod, rf3_chiral_off_baseline
util.cnc rf3_chiral_off_baseline
show cartoon, rf3_chiral_off_baseline

load ./5SOP_rf3_chiral_off_gw0p01.cif, rf3_chiral_off_gw_0_01
align rf3_chiral_off_gw_0_01, density_input
color good_goldenrod, rf3_chiral_off_gw_0_01
util.cnc rf3_chiral_off_gw_0_01
show cartoon, rf3_chiral_off_gw_0_01

# ---------------------------------------------------------------
# Groups for easy show/hide toggling (flat - one level only)
# ---------------------------------------------------------------
group boltz2_xray, boltz2_xray_baseline boltz2_xray_gw_0_1
group boltz2_md, boltz2_md_baseline boltz2_md_gw_0_1
group protenix, protenix_baseline protenix_gw_0_1
group protenix_mini, protenix_mini_baseline protenix_mini_gw_0_1
group protenix_tiny, protenix_tiny_baseline protenix_tiny_gw_0_1
group rf3_fixed, rf3_fixed_baseline rf3_fixed_gw_0_01
group rf3_chiral_off, rf3_chiral_off_baseline rf3_chiral_off_gw_0_01

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

# Per-residue display - resi 132 sticks on every object, hide resi 45-52
show sticks, resi 132
hide cartoon, resi 45-52

# Transparency - models see-through, density_input stays opaque
set cartoon_transparency, 0.3
set cartoon_transparency, 0.0, density_input
set stick_transparency, 0.3
set stick_transparency, 0.0, density_input

# ---------------------------------------------------------------
# VIEWING ANGLE - fixed view captured via get_view
# ---------------------------------------------------------------
set_view (\
     0.326664031,    0.585407972,   -0.742011845,\
    -0.885843515,   -0.084069811,   -0.456308991,\
    -0.329508662,    0.806362271,    0.491118968,\
    -0.000075484,   -0.000238800,  -44.233943939,\
    23.471965790,   21.350997925,   12.345628738,\
    31.563886642,   56.943958282,   20.000000000 )

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
enable protenix_mini
enable protenix_tiny
enable rf3_fixed
enable rf3_chiral_off

# --- Boltz2 X-ray ---
enable boltz2_xray_baseline
ray 2400, 1800
png 5SOP_boltz2_xray_baseline.png, dpi=300
disable boltz2_xray_baseline
enable boltz2_xray_gw_0_1
ray 2400, 1800
png 5SOP_boltz2_xray_guided.png, dpi=300
disable boltz2_xray_gw_0_1

# --- Boltz2 MD ---
enable boltz2_md_baseline
ray 2400, 1800
png 5SOP_boltz2_md_baseline.png, dpi=300
disable boltz2_md_baseline
enable boltz2_md_gw_0_1
ray 2400, 1800
png 5SOP_boltz2_md_guided.png, dpi=300
disable boltz2_md_gw_0_1

# --- Protenix ---
enable protenix_baseline
ray 2400, 1800
png 5SOP_protenix_baseline.png, dpi=300
disable protenix_baseline
enable protenix_gw_0_1
ray 2400, 1800
png 5SOP_protenix_guided.png, dpi=300
disable protenix_gw_0_1

# --- Protenix mini ---
enable protenix_mini_baseline
ray 2400, 1800
png 5SOP_protenix_mini_baseline.png, dpi=300
disable protenix_mini_baseline
enable protenix_mini_gw_0_1
ray 2400, 1800
png 5SOP_protenix_mini_guided.png, dpi=300
disable protenix_mini_gw_0_1

# --- Protenix tiny ---
enable protenix_tiny_baseline
ray 2400, 1800
png 5SOP_protenix_tiny_baseline.png, dpi=300
disable protenix_tiny_baseline
enable protenix_tiny_gw_0_1
ray 2400, 1800
png 5SOP_protenix_tiny_guided.png, dpi=300
disable protenix_tiny_gw_0_1

# --- RF3 fixed ---
enable rf3_fixed_baseline
ray 2400, 1800
png 5SOP_rf3_fixed_baseline.png, dpi=300
disable rf3_fixed_baseline
enable rf3_fixed_gw_0_01
ray 2400, 1800
png 5SOP_rf3_fixed_guided.png, dpi=300
disable rf3_fixed_gw_0_01

# --- RF3 chiral-off ---
enable rf3_chiral_off_baseline
ray 2400, 1800
png 5SOP_rf3_chiral_off_baseline.png, dpi=300
disable rf3_chiral_off_baseline
enable rf3_chiral_off_gw_0_01
ray 2400, 1800
png 5SOP_rf3_chiral_off_guided.png, dpi=300
disable rf3_chiral_off_gw_0_01

# Re-enable everything for interactive inspection
enable all
