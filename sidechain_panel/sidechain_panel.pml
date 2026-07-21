# Cross-model side-chain guidance panel - one protein per model
# Each model is shown on the residue that best demonstrates its side-chain
# guidance, baseline (gw0.0) vs guided (per-model optimal gw)
#   boltz2 X-ray : 3DT2 resi 11  PHE   gw0.0 -> gw0.1   good_darkblue
#   boltz2 MD    : 2QCV resi 261 LYS   gw0.0 -> gw0.1   good_teal
#   rf3          : 3L4P resi 18  ASP   gw0.0 -> gw0.01  good_navaho
#   protenix     : 3T94 resi 229 THR   gw0.0 -> gw0.1   good_pink
# Unlike the single-protein panels each model has its OWN protein, density
# input and coordinate frame, so they cannot share one view - every protein is
# framed and rendered independently over its own gray density reference
# Each model object is an ens8 stack of 8 states - baseline and guided share one
# color and are distinguished by toggling show/hide before each ray

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

# ===============================================================
# Boltz2 X-ray   3DT2 resi 11 PHE   (good_darkblue)
# ===============================================================
load_altfix ./3DT2_density_input.cif, boltz2_xray_density
color gray70, boltz2_xray_density
util.cnc boltz2_xray_density
show cartoon, boltz2_xray_density

load ./3DT2_boltz2_xray_gw0p0.cif, boltz2_xray_baseline
color good_darkblue, boltz2_xray_baseline
util.cnc boltz2_xray_baseline
show cartoon, boltz2_xray_baseline

load ./3DT2_boltz2_xray_gw0p1.cif, boltz2_xray_guided
color good_darkblue, boltz2_xray_guided
util.cnc boltz2_xray_guided
show cartoon, boltz2_xray_guided

# ===============================================================
# Boltz2 MD   2QCV resi 261 LYS   (good_teal)
# ===============================================================
load_altfix ./2QCV_density_input.cif, boltz2_md_density
color gray70, boltz2_md_density
util.cnc boltz2_md_density
show cartoon, boltz2_md_density

load ./2QCV_boltz2_md_gw0p0.cif, boltz2_md_baseline
color good_teal, boltz2_md_baseline
util.cnc boltz2_md_baseline
show cartoon, boltz2_md_baseline

load ./2QCV_boltz2_md_gw0p1.cif, boltz2_md_guided
color good_teal, boltz2_md_guided
util.cnc boltz2_md_guided
show cartoon, boltz2_md_guided

# ===============================================================
# RF3 fixed   3L4P resi 18 ASP   (good_navaho)
# ===============================================================
load_altfix ./3L4P_density_input.cif, rf3_fixed_density
color gray70, rf3_fixed_density
util.cnc rf3_fixed_density
show cartoon, rf3_fixed_density

load ./3L4P_rf3_fixed_gw0p0.cif, rf3_fixed_baseline
color good_navaho, rf3_fixed_baseline
util.cnc rf3_fixed_baseline
show cartoon, rf3_fixed_baseline

load ./3L4P_rf3_fixed_gw0p01.cif, rf3_fixed_guided
color good_navaho, rf3_fixed_guided
util.cnc rf3_fixed_guided
show cartoon, rf3_fixed_guided

# ===============================================================
# Protenix   3T94 resi 229 THR   (good_pink)
# ===============================================================
load_altfix ./3T94_density_input.cif, protenix_density
color gray70, protenix_density
util.cnc protenix_density
show cartoon, protenix_density

load ./3T94_protenix_gw0p0.cif, protenix_baseline
color good_pink, protenix_baseline
util.cnc protenix_baseline
show cartoon, protenix_baseline

load ./3T94_protenix_gw0p1.cif, protenix_guided
color good_pink, protenix_guided
util.cnc protenix_guided
show cartoon, protenix_guided

# ---------------------------------------------------------------
# Groups for easy show/hide toggling (flat - one level only)
# Density objects are kept ungrouped so they toggle independently
# ---------------------------------------------------------------
group boltz2_xray, boltz2_xray_baseline boltz2_xray_guided
group boltz2_md, boltz2_md_baseline boltz2_md_guided
group rf3_fixed, rf3_fixed_baseline rf3_fixed_guided
group protenix, protenix_baseline protenix_guided

# ---------------------------------------------------------------
# Display settings
# ---------------------------------------------------------------
# Hide chain-break dashes
set cartoon_gap_cutoff, 0
# Render all model cartoons as loop - removes helix and sheet secondary structure
# Density objects keep their normal cartoon so the reference reads as the crystal
cartoon loop, not (boltz2_xray_density or boltz2_md_density or rf3_fixed_density or protenix_density)
# Background and depth - override .pymolrc as needed
bg_color white
set ray_shadow, 0
set ambient_occlusion_mode, 0
# Cartoon outlines (BLOPIG method) - normal color plus an edge line during ray
# ray_trace_mode 1 = color + outline. Higher ray_trace_gain = thicker/darker edge.
# At 2400x1800 the near-default 0.1 renders a sub-pixel hairline, so push it up.
# Outline color is inherited from .pymolrc (gray30) - antialias 2 also from there.
# Transparent objects receive NO outline in mode 1, so only the opaque gray density
# references are edged - the see-through model cartoons stay unoutlined by design.
set ray_trace_mode, 1
set ray_trace_gain, 1.5
# Show all 8 ensemble states at once
set state, 0, all

# Per-residue side-chain sticks - scoped per protein so each residue number only
# lights up on its own object set (resi 11 etc. also exist on other proteins)
show sticks, (boltz2_xray_density or boltz2_xray_baseline or boltz2_xray_guided) and chain A and resi 11
show sticks, (boltz2_md_density or boltz2_md_baseline or boltz2_md_guided) and chain A and resi 261
show sticks, (rf3_fixed_density or rf3_fixed_baseline or rf3_fixed_guided) and chain A and resi 18
show sticks, (protenix_density or protenix_baseline or protenix_guided) and chain A and resi 229

# Transparency - models see-through, density references stay opaque
set cartoon_transparency, 0.3
set cartoon_transparency, 0.0, boltz2_xray_density
set cartoon_transparency, 0.0, boltz2_md_density
set cartoon_transparency, 0.0, rf3_fixed_density
set cartoon_transparency, 0.0, protenix_density
set stick_transparency, 0.3
set stick_transparency, 0.0, boltz2_xray_density
set stick_transparency, 0.0, boltz2_md_density
set stick_transparency, 0.0, rf3_fixed_density
set stick_transparency, 0.0, protenix_density

# ---------------------------------------------------------------
# INTERACTIVE INSPECTION  (ray / png output disabled)
# Use this to judge each protein's framing and to spot any cartoon or stick that
# obstructs the residue before committing to full renders. Each block switches
# the scene to one protein - its gray density plus both ensemble stacks - and
# frames its residue. The four proteins live in different coordinate frames so
# only one can be viewed at a time. Sourcing the whole script lands on the last
# block (protenix) - run an earlier block to jump to that protein instead.
#
# To pin a framing - fine-tune the angle in the GUI, run get_view, and paste the
# 18-number set_view ( ... ) matrix just below that block's zoom line.
#
# To RESTORE RENDERING - in each block uncomment the lines that start with
#   ray / png / disable
# (one ray and png per condition, with the disable that hides it before the next
# condition). That returns to one-condition-at-a-time output to the PNGs below.
# ---------------------------------------------------------------

# --- Boltz2 X-ray   3DT2 resi 11   (baseline gw0.0 -> guided gw0.1) ---
disable all
enable boltz2_xray_density
enable boltz2_xray
enable boltz2_xray_baseline

set_view (\
     0.892827809,   -0.089483961,    0.441403270,\
     0.070224486,   -0.940416694,   -0.332685739,\
     0.444876641,    0.328026205,   -0.833355844,\
    -0.000682462,   -0.000059779,  -30.969341278,\
    67.636505127,   35.955684662,   19.360940933,\
     6.262629032,   55.589450836,   20.000000000 )

ray 2400, 1800
png 3DT2_boltz2_xray_baseline.png, dpi=300
disable boltz2_xray_baseline
enable boltz2_xray_guided
ray 2400, 1800
png 3DT2_boltz2_xray_guided.png, dpi=300
disable boltz2_xray_guided

# --- Boltz2 MD   2QCV resi 261   (baseline gw0.0 -> guided gw0.1) ---
disable all
enable boltz2_md_density
enable boltz2_md

set_view (\
     0.814386606,   -0.526949644,   -0.243100658,\
     0.093570389,    0.532657325,   -0.841142356,\
     0.572727859,    0.662270069,    0.483095765,\
    -0.000016080,    0.000041772,  -35.712059021,\
    12.011948586,   40.138881683,   43.601547241,\
    23.424983978,   48.002113342,   20.000000000 )

enable boltz2_md_baseline
ray 2400, 1800
png 2QCV_boltz2_md_baseline.png, dpi=300
disable boltz2_md_baseline
enable boltz2_md_guided
ray 2400, 1800
png 2QCV_boltz2_md_guided.png, dpi=300
disable boltz2_md_guided

# --- RF3 fixed   3L4P resi 18   (baseline gw0.0 -> guided gw0.01) ---
disable all
enable rf3_fixed_density
enable rf3_fixed

set_view (\
    -0.606252730,   -0.578937173,   -0.545242906,\
    -0.719259083,    0.106641270,    0.686508656,\
    -0.339299023,    0.808368623,   -0.481055588,\
     0.000060770,   -0.000025185,  -23.480501175,\
    33.521980286,   67.163902283,   23.487360001,\
    20.245923996,   26.716207504,   20.000000000 )

enable rf3_fixed_baseline
ray 2400, 1800
png 3L4P_rf3_fixed_baseline.png, dpi=300
disable rf3_fixed_baseline
enable rf3_fixed_guided
ray 2400, 1800
png 3L4P_rf3_fixed_guided.png, dpi=300
disable rf3_fixed_guided

# --- Protenix   3T94 resi 229   (baseline gw0.0 -> guided gw0.1) ---
disable all
enable protenix_density
enable protenix

set_view (\
     0.845125973,    0.337502122,    0.414544702,\
    -0.070633903,    0.839190602,   -0.539229989,\
    -0.529875994,    0.426435441,    0.733064055,\
    -0.000076956,   -0.000254339,  -30.407516479,\
    17.343235016,   37.955944061,   37.073310852,\
    -3.775797129,   64.541015625,   20.000000000 )

enable protenix_baseline
ray 2400, 1800
png 3T94_protenix_baseline.png, dpi=300
disable protenix_baseline
enable protenix_guided
ray 2400, 1800
png 3T94_protenix_guided.png, dpi=300

