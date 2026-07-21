# 5I09 0.5occA_0.5occB panel - geometry decay across the guidance sweep
# 4 MAIN models only - boltz2 X-ray, boltz2 MD, protenix, rf3
# Each model shown at 4 guidance weights: baseline -> mid -> optimal -> highest
#   boltz2 / protenix : gw0.0  gw0.05   gw0.1   gw0.5
#   rf3               : gw0.0  gw0.005  gw0.01  gw0.1
# Goal: show geometry degrading as guidance weight climbs past the optimum
# Each model object is an ens8 stack (8 states), all 4 gw levels share one color
# Distinguish the gw levels by toggling show/hide before each ray

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
# Density input (reference) - stays in the map frame
# ---------------------------------------------------------------
load_altfix ./5I09_density_input.cif, density_input
color gray70, density_input
util.cnc density_input
show cartoon, density_input

# ---------------------------------------------------------------
# Boltz2 X-ray density   (good_darkblue)
# ---------------------------------------------------------------
load ./5I09_boltz2_xray_gw0p0.cif, boltz2_xray_baseline
align boltz2_xray_baseline, density_input
color good_darkblue, boltz2_xray_baseline
util.cnc boltz2_xray_baseline
show cartoon, boltz2_xray_baseline

load ./5I09_boltz2_xray_gw0p05.cif, boltz2_xray_mid
align boltz2_xray_mid, density_input
color good_darkblue, boltz2_xray_mid
util.cnc boltz2_xray_mid
show cartoon, boltz2_xray_mid

load ./5I09_boltz2_xray_gw0p1.cif, boltz2_xray_optimal
align boltz2_xray_optimal, density_input
color good_darkblue, boltz2_xray_optimal
util.cnc boltz2_xray_optimal
show cartoon, boltz2_xray_optimal

load ./5I09_boltz2_xray_gw0p5.cif, boltz2_xray_highest
align boltz2_xray_highest, density_input
color good_darkblue, boltz2_xray_highest
util.cnc boltz2_xray_highest
show cartoon, boltz2_xray_highest

# ---------------------------------------------------------------
# Boltz2 MD density   (good_teal)
# ---------------------------------------------------------------
load ./5I09_boltz2_md_gw0p0.cif, boltz2_md_baseline
align boltz2_md_baseline, density_input
color good_teal, boltz2_md_baseline
util.cnc boltz2_md_baseline
show cartoon, boltz2_md_baseline

load ./5I09_boltz2_md_gw0p05.cif, boltz2_md_mid
align boltz2_md_mid, density_input
color good_teal, boltz2_md_mid
util.cnc boltz2_md_mid
show cartoon, boltz2_md_mid

load ./5I09_boltz2_md_gw0p1.cif, boltz2_md_optimal
align boltz2_md_optimal, density_input
color good_teal, boltz2_md_optimal
util.cnc boltz2_md_optimal
show cartoon, boltz2_md_optimal

load ./5I09_boltz2_md_gw0p5.cif, boltz2_md_highest
align boltz2_md_highest, density_input
color good_teal, boltz2_md_highest
util.cnc boltz2_md_highest
show cartoon, boltz2_md_highest

# ---------------------------------------------------------------
# Protenix   (good_pink)
# ---------------------------------------------------------------
load ./5I09_protenix_gw0p0.cif, protenix_baseline
align protenix_baseline, density_input
color good_pink, protenix_baseline
util.cnc protenix_baseline
show cartoon, protenix_baseline

load ./5I09_protenix_gw0p05.cif, protenix_mid
align protenix_mid, density_input
color good_pink, protenix_mid
util.cnc protenix_mid
show cartoon, protenix_mid

load ./5I09_protenix_gw0p1.cif, protenix_optimal
align protenix_optimal, density_input
color good_pink, protenix_optimal
util.cnc protenix_optimal
show cartoon, protenix_optimal

load ./5I09_protenix_gw0p5.cif, protenix_highest
align protenix_highest, density_input
color good_pink, protenix_highest
util.cnc protenix_highest
show cartoon, protenix_highest

# ---------------------------------------------------------------
# RF3 fixed   (good_navaho)
# ---------------------------------------------------------------
load ./5I09_rf3_fixed_gw0p0.cif, rf3_fixed_baseline
align rf3_fixed_baseline, density_input
color good_navaho, rf3_fixed_baseline
util.cnc rf3_fixed_baseline
show cartoon, rf3_fixed_baseline

load ./5I09_rf3_fixed_gw0p005.cif, rf3_fixed_mid
align rf3_fixed_mid, density_input
color good_navaho, rf3_fixed_mid
util.cnc rf3_fixed_mid
show cartoon, rf3_fixed_mid

load ./5I09_rf3_fixed_gw0p01.cif, rf3_fixed_optimal
align rf3_fixed_optimal, density_input
color good_navaho, rf3_fixed_optimal
util.cnc rf3_fixed_optimal
show cartoon, rf3_fixed_optimal

load ./5I09_rf3_fixed_gw0p1.cif, rf3_fixed_highest
align rf3_fixed_highest, density_input
color good_navaho, rf3_fixed_highest
util.cnc rf3_fixed_highest
show cartoon, rf3_fixed_highest

# ---------------------------------------------------------------
# Groups for easy show/hide toggling (flat - one level only)
# ---------------------------------------------------------------
group boltz2_xray, boltz2_xray_baseline boltz2_xray_mid boltz2_xray_optimal boltz2_xray_highest
group boltz2_md, boltz2_md_baseline boltz2_md_mid boltz2_md_optimal boltz2_md_highest
group protenix, protenix_baseline protenix_mid protenix_optimal protenix_highest
group rf3_fixed, rf3_fixed_baseline rf3_fixed_mid rf3_fixed_optimal rf3_fixed_highest

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
# Outline color is inherited from .pymolrc (gray30), antialias 2 also from there.
# Note: transparent objects receive NO outline in mode 1, so only opaque objects
# are edged (here that is density_input, see the transparency block below).
set ray_trace_mode, 1
set ray_trace_gain, 1.5
# Show all 8 ensemble states at once
set state, 0, all

# Per-residue focus (OPTIONAL) - 5I09 has no single dominant altloc region, so
# the whole structure is shown and geometry decay reads globally as the ens8
# stack spreads apart. To spotlight a region instead, set a view below and
# uncomment/edit lines like these (the largest 5I09 backbone altloc is 68-73):
show sticks, resi 68-74
hide cartoon, resi 68-74
hide everything, resi 43-54 or resi 280-283

# Transparency - models see-through, density_input stays opaque
# Bump a model to 0.0 to make it opaque and pick up a ray-trace outline.
set cartoon_transparency, 0.3
set cartoon_transparency, 0.0, density_input
set stick_transparency, 0.3
set stick_transparency, 0.0, density_input

# ---------------------------------------------------------------
# VIEWING ANGLE - PLACEHOLDER (replace with a captured view)
# Open this script in PyMOL, rotate/zoom to a good angle, then run
#   get_view
# and paste the 18-number matrix here as a set_view ( ... ) block,
# removing the orient line below.
# ---------------------------------------------------------------
set_view (\
    -0.515102267,    0.816999674,    0.259168863,\
     0.693674922,    0.219735637,    0.685947299,\
     0.503468513,    0.533110201,   -0.679920197,\
    -0.000609387,   -0.001411408,  -32.582977295,\
    22.066518784,   37.054958344,   45.109432220,\
    21.320953369,   43.895206451,   20.000000000 )


# ---------------------------------------------------------------
# RAY TRACE / EXPORT
# Renders one model gw-level at a time over the gray reference
# Density input stays visible throughout
# Comment out any model block you do not need - each ray is full-res
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

# --- Boltz2 X-ray  (baseline gw0.0 -> mid gw0.05 -> optimal gw0.1 -> highest gw0.5) ---
enable boltz2_xray_baseline
ray 2400, 1800
png 5I09_boltz2_xray_baseline.png, dpi=300
disable boltz2_xray_baseline
enable boltz2_xray_mid
ray 2400, 1800
png 5I09_boltz2_xray_mid.png, dpi=300
disable boltz2_xray_mid
enable boltz2_xray_optimal
ray 2400, 1800
png 5I09_boltz2_xray_optimal.png, dpi=300
disable boltz2_xray_optimal
enable boltz2_xray_highest
ray 2400, 1800
png 5I09_boltz2_xray_highest.png, dpi=300
disable boltz2_xray_highest

# --- Boltz2 MD  (baseline gw0.0 -> mid gw0.05 -> optimal gw0.1 -> highest gw0.5) ---
enable boltz2_md_baseline
ray 2400, 1800
png 5I09_boltz2_md_baseline.png, dpi=300
disable boltz2_md_baseline
enable boltz2_md_mid
ray 2400, 1800
png 5I09_boltz2_md_mid.png, dpi=300
disable boltz2_md_mid
enable boltz2_md_optimal
ray 2400, 1800
png 5I09_boltz2_md_optimal.png, dpi=300
disable boltz2_md_optimal
enable boltz2_md_highest
ray 2400, 1800
png 5I09_boltz2_md_highest.png, dpi=300
disable boltz2_md_highest

# --- Protenix  (baseline gw0.0 -> mid gw0.05 -> optimal gw0.1 -> highest gw0.5) ---
enable protenix_baseline
ray 2400, 1800
png 5I09_protenix_baseline.png, dpi=300
disable protenix_baseline
enable protenix_mid
ray 2400, 1800
png 5I09_protenix_mid.png, dpi=300
disable protenix_mid
enable protenix_optimal
ray 2400, 1800
png 5I09_protenix_optimal.png, dpi=300
disable protenix_optimal
enable protenix_highest
ray 2400, 1800
png 5I09_protenix_highest.png, dpi=300
disable protenix_highest

# --- RF3 fixed  (baseline gw0.0 -> mid gw0.005 -> optimal gw0.01 -> highest gw0.1) ---
enable rf3_fixed_baseline
ray 2400, 1800
png 5I09_rf3_fixed_baseline.png, dpi=300
disable rf3_fixed_baseline
enable rf3_fixed_mid
ray 2400, 1800
png 5I09_rf3_fixed_mid.png, dpi=300
disable rf3_fixed_mid
enable rf3_fixed_optimal
ray 2400, 1800
png 5I09_rf3_fixed_optimal.png, dpi=300
disable rf3_fixed_optimal
enable rf3_fixed_highest
ray 2400, 1800
png 5I09_rf3_fixed_highest.png, dpi=300
disable rf3_fixed_highest

# Re-enable everything for interactive inspection
enable all
