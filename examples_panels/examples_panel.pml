# Examples panel - 6B8X, 5MHX, 5SOP baseline-vs-guided ensembles vs experimental density
# Single self-contained port of baseline/view.pml. All inputs copied locally with
# flat names by copy_examples_cifs.sh (run from this folder:  pymol examples_panel.pml).
#
# 4 MAIN models per structure - boltz2 X-ray, boltz2 MD, protenix, rf3 - each shown
# baseline (gw0.0) vs guided (optimal gw):
#   boltz2 / protenix : gw0.0  gw0.1
#   rf3               : gw0.0  gw0.01
# Each model object is a refined-patched ens8 stack (8 states). Baseline and guided
# share one color - distinguish them by toggling show/hide before each ray.
# A gray density_input (keeps crystallographic alt A + alt B) is the reference, and a
# carved isomesh of the experimental 1.00 A map (slate) shows the target density.

# ---------------------------------------------------------------
# Custom colors - mirror ~/.pymolrc so this script is self-contained
# ---------------------------------------------------------------
set_color good_darkblue, [75, 95, 170]
set_color good_teal,     [0.310, 0.725, 0.686]
set_color good_navaho,   [255, 224, 172]
set_color good_pink,     [255, 172, 183]

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
# 6B8X
# ===============================================================
load_altfix ./6B8X_gt.cif, 6B8X_gt
color gray70, 6B8X_gt
util.cnc 6B8X_gt
show cartoon, 6B8X_gt
load ./6B8X_map.ccp4, 6B8X_map

# boltz2 X-ray (good_darkblue)
load ./6B8X_boltz2XRAY_gw0p0.cif, 6B8X_boltz2XRAY_baseline
align 6B8X_boltz2XRAY_baseline, 6B8X_gt
color good_darkblue, 6B8X_boltz2XRAY_baseline
util.cnc 6B8X_boltz2XRAY_baseline
show cartoon, 6B8X_boltz2XRAY_baseline
load ./6B8X_boltz2XRAY_gw0p1.cif, 6B8X_boltz2XRAY_guided
align 6B8X_boltz2XRAY_guided, 6B8X_gt
color good_darkblue, 6B8X_boltz2XRAY_guided
util.cnc 6B8X_boltz2XRAY_guided
show cartoon, 6B8X_boltz2XRAY_guided

# boltz2 MD (good_teal)
load ./6B8X_boltz2MD_gw0p0.cif, 6B8X_boltz2MD_baseline
align 6B8X_boltz2MD_baseline, 6B8X_gt
color good_teal, 6B8X_boltz2MD_baseline
util.cnc 6B8X_boltz2MD_baseline
show cartoon, 6B8X_boltz2MD_baseline
load ./6B8X_boltz2MD_gw0p1.cif, 6B8X_boltz2MD_guided
align 6B8X_boltz2MD_guided, 6B8X_gt
color good_teal, 6B8X_boltz2MD_guided
util.cnc 6B8X_boltz2MD_guided
show cartoon, 6B8X_boltz2MD_guided

# protenix (good_pink)
load ./6B8X_protenix_gw0p0.cif, 6B8X_protenix_baseline
align 6B8X_protenix_baseline, 6B8X_gt
color good_pink, 6B8X_protenix_baseline
util.cnc 6B8X_protenix_baseline
show cartoon, 6B8X_protenix_baseline
load ./6B8X_protenix_gw0p1.cif, 6B8X_protenix_guided
align 6B8X_protenix_guided, 6B8X_gt
color good_pink, 6B8X_protenix_guided
util.cnc 6B8X_protenix_guided
show cartoon, 6B8X_protenix_guided

# rf3 (good_navaho)
load ./6B8X_rf3_gw0p0.cif, 6B8X_rf3_baseline
align 6B8X_rf3_baseline, 6B8X_gt
color good_navaho, 6B8X_rf3_baseline
util.cnc 6B8X_rf3_baseline
show cartoon, 6B8X_rf3_baseline
load ./6B8X_rf3_gw0p01.cif, 6B8X_rf3_guided
align 6B8X_rf3_guided, 6B8X_gt
color good_navaho, 6B8X_rf3_guided
util.cnc 6B8X_rf3_guided
show cartoon, 6B8X_rf3_guided

# ===============================================================
# 5MHX
# ===============================================================
load_altfix ./5MHX_gt.cif, 5MHX_gt
color gray70, 5MHX_gt
util.cnc 5MHX_gt
show cartoon, 5MHX_gt
load ./5MHX_map.ccp4, 5MHX_map

# boltz2 X-ray (good_darkblue)
load ./5MHX_boltz2XRAY_gw0p0.cif, 5MHX_boltz2XRAY_baseline
align 5MHX_boltz2XRAY_baseline, 5MHX_gt
color good_darkblue, 5MHX_boltz2XRAY_baseline
util.cnc 5MHX_boltz2XRAY_baseline
show cartoon, 5MHX_boltz2XRAY_baseline
load ./5MHX_boltz2XRAY_gw0p1.cif, 5MHX_boltz2XRAY_guided
align 5MHX_boltz2XRAY_guided, 5MHX_gt
color good_darkblue, 5MHX_boltz2XRAY_guided
util.cnc 5MHX_boltz2XRAY_guided
show cartoon, 5MHX_boltz2XRAY_guided

# boltz2 MD (good_teal)
load ./5MHX_boltz2MD_gw0p0.cif, 5MHX_boltz2MD_baseline
align 5MHX_boltz2MD_baseline, 5MHX_gt
color good_teal, 5MHX_boltz2MD_baseline
util.cnc 5MHX_boltz2MD_baseline
show cartoon, 5MHX_boltz2MD_baseline
load ./5MHX_boltz2MD_gw0p1.cif, 5MHX_boltz2MD_guided
align 5MHX_boltz2MD_guided, 5MHX_gt
color good_teal, 5MHX_boltz2MD_guided
util.cnc 5MHX_boltz2MD_guided
show cartoon, 5MHX_boltz2MD_guided

# protenix (good_pink)
load ./5MHX_protenix_gw0p0.cif, 5MHX_protenix_baseline
align 5MHX_protenix_baseline, 5MHX_gt
color good_pink, 5MHX_protenix_baseline
util.cnc 5MHX_protenix_baseline
show cartoon, 5MHX_protenix_baseline
load ./5MHX_protenix_gw0p1.cif, 5MHX_protenix_guided
align 5MHX_protenix_guided, 5MHX_gt
color good_pink, 5MHX_protenix_guided
util.cnc 5MHX_protenix_guided
show cartoon, 5MHX_protenix_guided

# rf3 (good_navaho)
load ./5MHX_rf3_gw0p0.cif, 5MHX_rf3_baseline
align 5MHX_rf3_baseline, 5MHX_gt
color good_navaho, 5MHX_rf3_baseline
util.cnc 5MHX_rf3_baseline
show cartoon, 5MHX_rf3_baseline
load ./5MHX_rf3_gw0p01.cif, 5MHX_rf3_guided
align 5MHX_rf3_guided, 5MHX_gt
color good_navaho, 5MHX_rf3_guided
util.cnc 5MHX_rf3_guided
show cartoon, 5MHX_rf3_guided

# ===============================================================
# 5SOP   (density_input has a spurious alt C - drop it)
# ===============================================================
load_altfix ./5SOP_gt.cif, 5SOP_gt
remove 5SOP_gt and alt C
color gray70, 5SOP_gt
util.cnc 5SOP_gt
show cartoon, 5SOP_gt
load ./5SOP_map.ccp4, 5SOP_map

# boltz2 X-ray (good_darkblue)
load ./5SOP_boltz2XRAY_gw0p0.cif, 5SOP_boltz2XRAY_baseline
align 5SOP_boltz2XRAY_baseline, 5SOP_gt
color good_darkblue, 5SOP_boltz2XRAY_baseline
util.cnc 5SOP_boltz2XRAY_baseline
show cartoon, 5SOP_boltz2XRAY_baseline
load ./5SOP_boltz2XRAY_gw0p1.cif, 5SOP_boltz2XRAY_guided
align 5SOP_boltz2XRAY_guided, 5SOP_gt
color good_darkblue, 5SOP_boltz2XRAY_guided
util.cnc 5SOP_boltz2XRAY_guided
show cartoon, 5SOP_boltz2XRAY_guided

# boltz2 MD (good_teal)
load ./5SOP_boltz2MD_gw0p0.cif, 5SOP_boltz2MD_baseline
align 5SOP_boltz2MD_baseline, 5SOP_gt
color good_teal, 5SOP_boltz2MD_baseline
util.cnc 5SOP_boltz2MD_baseline
show cartoon, 5SOP_boltz2MD_baseline
load ./5SOP_boltz2MD_gw0p1.cif, 5SOP_boltz2MD_guided
align 5SOP_boltz2MD_guided, 5SOP_gt
color good_teal, 5SOP_boltz2MD_guided
util.cnc 5SOP_boltz2MD_guided
show cartoon, 5SOP_boltz2MD_guided

# protenix (good_pink)
load ./5SOP_protenix_gw0p0.cif, 5SOP_protenix_baseline
align 5SOP_protenix_baseline, 5SOP_gt
color good_pink, 5SOP_protenix_baseline
util.cnc 5SOP_protenix_baseline
show cartoon, 5SOP_protenix_baseline
load ./5SOP_protenix_gw0p1.cif, 5SOP_protenix_guided
align 5SOP_protenix_guided, 5SOP_gt
color good_pink, 5SOP_protenix_guided
util.cnc 5SOP_protenix_guided
show cartoon, 5SOP_protenix_guided

# rf3 (good_navaho)
load ./5SOP_rf3_gw0p0.cif, 5SOP_rf3_baseline
align 5SOP_rf3_baseline, 5SOP_gt
color good_navaho, 5SOP_rf3_baseline
util.cnc 5SOP_rf3_baseline
show cartoon, 5SOP_rf3_baseline
load ./5SOP_rf3_gw0p01.cif, 5SOP_rf3_guided
align 5SOP_rf3_guided, 5SOP_gt
color good_navaho, 5SOP_rf3_guided
util.cnc 5SOP_rf3_guided
show cartoon, 5SOP_rf3_guided

# Maps are only the isomesh source - keep them off
disable *_map

# ===============================================================
# Groups (flat, one level) - per model: baseline + guided
# ===============================================================
group 6B8X_boltz2XRAY, 6B8X_boltz2XRAY_baseline 6B8X_boltz2XRAY_guided
group 6B8X_boltz2MD,   6B8X_boltz2MD_baseline 6B8X_boltz2MD_guided
group 6B8X_protenix,   6B8X_protenix_baseline 6B8X_protenix_guided
group 6B8X_rf3,        6B8X_rf3_baseline 6B8X_rf3_guided

group 5MHX_boltz2XRAY, 5MHX_boltz2XRAY_baseline 5MHX_boltz2XRAY_guided
group 5MHX_boltz2MD,   5MHX_boltz2MD_baseline 5MHX_boltz2MD_guided
group 5MHX_protenix,   5MHX_protenix_baseline 5MHX_protenix_guided
group 5MHX_rf3,        5MHX_rf3_baseline 5MHX_rf3_guided

group 5SOP_boltz2XRAY, 5SOP_boltz2XRAY_baseline 5SOP_boltz2XRAY_guided
group 5SOP_boltz2MD,   5SOP_boltz2MD_baseline 5SOP_boltz2MD_guided
group 5SOP_protenix,   5SOP_protenix_baseline 5SOP_protenix_guided
group 5SOP_rf3,        5SOP_rf3_baseline 5SOP_rf3_guided

# ===============================================================
# Display settings
# ===============================================================
set cartoon_transparency, 0.3
set cartoon_transparency, 0.0, *_gt
set stick_transparency, 0.3
set stick_transparency, 0.0, *_gt
bg_color white
set ray_shadow, 0
set ambient_occlusion_mode, 0
# Black outlines on opaque GT via ray_trace_mode 1 (BLOPIG method).
# Transparent model cartoons receive NO outline in this mode (only opaque GT).
set ray_trace_mode, 1
set ray_trace_gain, 0.1
set cartoon_gap_cutoff, 0
# Show all 8 ensemble states at once
set state, 0

# ===============================================================
# Density isomesh (carved from the experimental 1.00 A map around each altloc)
# ===============================================================
isomesh 6B8X_isomesh, 6B8X_map, 1.0, 6B8X_gt and resi 175-184, carve=1.7
isomesh 5MHX_isomesh, 5MHX_map, 1.0, 5MHX_gt and resi 158-165, carve=1.8
isomesh 5SOP_isomesh, 5SOP_map, 1.0, 5SOP_gt and resi 127-136, carve=1.8
hide everything, resi 47-49 and 5SOP*
color slate, *_isomesh
set transparency, 0.4, *_isomesh

# ===============================================================
# Stick selections (one per structure, span GT + all models)
# Refined-patched CIFs use uniform numbering, so 5SOP needs no shifted select.
# ===============================================================
select select_6b8x, resi 181-183 and 6B8X*
select select_5mhx, resi 158-165 and 5MHX*
select select_5sop, resi 132 and 5SOP*

# ===============================================================
# RAY TRACE / EXPORT
# Antialias OFF - required for clean ray_trace_mode 1 outlines.
# Per structure: pin the view, render the density blob, then each model's
# baseline -> guided over the gray reference. Density input stays visible throughout.
# ===============================================================
set antialias, 0

# ---------------------------------------------------------------
# 6B8X shots
# ---------------------------------------------------------------
disable
enable 6B8X_gt
# Un-mask the model group containers so enabling a leaf below renders it
enable 6B8X_boltz2XRAY
enable 6B8X_boltz2MD
enable 6B8X_protenix
enable 6B8X_rf3
show sticks, select_6b8x
# Hide GT broken-chain stretches
hide everything, 6B8X_gt and (resi 108-114 or resi 116-119 or resi 215 or resi 217)

set_view (\
     0.767313063,    0.062839039,   -0.638187468,\
    -0.520306230,    0.642734170,   -0.562293351,\
     0.374851197,    0.763508201,    0.525872827,\
     0.000075649,    0.000056289,  -46.990047455,\
    31.870944977,   39.426567078,   41.271205902,\
    33.964267731,   60.015323639,  -20.000000000 )

# GT + density mesh
enable 6B8X_isomesh
ray 1680, 1262
png ./6B8X_density.png, dpi=300
disable 6B8X_isomesh

# boltz2 X-ray (baseline gw0.0 -> guided gw0.1)
enable 6B8X_boltz2XRAY_baseline
ray 1680, 1262
png ./6B8X_boltz2XRAY_baseline.png, dpi=300
disable 6B8X_boltz2XRAY_baseline
enable 6B8X_boltz2XRAY_guided
ray 1680, 1262
png ./6B8X_boltz2XRAY_guided.png, dpi=300
disable 6B8X_boltz2XRAY_guided

# boltz2 MD (baseline gw0.0 -> guided gw0.1)
enable 6B8X_boltz2MD_baseline
ray 1680, 1262
png ./6B8X_boltz2MD_baseline.png, dpi=300
disable 6B8X_boltz2MD_baseline
enable 6B8X_boltz2MD_guided
ray 1680, 1262
png ./6B8X_boltz2MD_guided.png, dpi=300
disable 6B8X_boltz2MD_guided

# protenix (baseline gw0.0 -> guided gw0.1)
enable 6B8X_protenix_baseline
ray 1680, 1262
png ./6B8X_protenix_baseline.png, dpi=300
disable 6B8X_protenix_baseline
enable 6B8X_protenix_guided
ray 1680, 1262
png ./6B8X_protenix_guided.png, dpi=300
disable 6B8X_protenix_guided

# rf3 (baseline gw0.0 -> guided gw0.01)
enable 6B8X_rf3_baseline
ray 1680, 1262
png ./6B8X_rf3_baseline.png, dpi=300
disable 6B8X_rf3_baseline
enable 6B8X_rf3_guided
ray 1680, 1262
png ./6B8X_rf3_guided.png, dpi=300
disable 6B8X_rf3_guided

# ---------------------------------------------------------------
# 5MHX shots
# ---------------------------------------------------------------
disable
enable 5MHX_gt
enable 5MHX_boltz2XRAY
enable 5MHX_boltz2MD
enable 5MHX_protenix
enable 5MHX_rf3
show sticks, select_5mhx

set_view (\
     0.666493356,   -0.241294548,   -0.705381513,\
    -0.070008390,    0.921731591,   -0.381453305,\
     0.742213428,    0.303619474,    0.597434521,\
    -0.000160942,   -0.000145167,  -36.099937439,\
    21.044620514,   47.041435242,   43.043750763,\
    19.171716690,   53.024841309,   20.000000000 )

# GT + density mesh
enable 5MHX_isomesh
ray 1680, 1262
png ./5MHX_density.png, dpi=300
disable 5MHX_isomesh

# boltz2 X-ray
enable 5MHX_boltz2XRAY_baseline
ray 1680, 1262
png ./5MHX_boltz2XRAY_baseline.png, dpi=300
disable 5MHX_boltz2XRAY_baseline
enable 5MHX_boltz2XRAY_guided
ray 1680, 1262
png ./5MHX_boltz2XRAY_guided.png, dpi=300
disable 5MHX_boltz2XRAY_guided

# boltz2 MD
enable 5MHX_boltz2MD_baseline
ray 1680, 1262
png ./5MHX_boltz2MD_baseline.png, dpi=300
disable 5MHX_boltz2MD_baseline
enable 5MHX_boltz2MD_guided
ray 1680, 1262
png ./5MHX_boltz2MD_guided.png, dpi=300
disable 5MHX_boltz2MD_guided

# protenix
enable 5MHX_protenix_baseline
ray 1680, 1262
png ./5MHX_protenix_baseline.png, dpi=300
disable 5MHX_protenix_baseline
enable 5MHX_protenix_guided
ray 1680, 1262
png ./5MHX_protenix_guided.png, dpi=300
disable 5MHX_protenix_guided

# rf3
enable 5MHX_rf3_baseline
ray 1680, 1262
png ./5MHX_rf3_baseline.png, dpi=300
disable 5MHX_rf3_baseline
enable 5MHX_rf3_guided
ray 1680, 1262
png ./5MHX_rf3_guided.png, dpi=300
disable 5MHX_rf3_guided

# ---------------------------------------------------------------
# 5SOP shots
# ---------------------------------------------------------------
disable
enable 5SOP_gt
enable 5SOP_boltz2XRAY
enable 5SOP_boltz2MD
enable 5SOP_protenix
enable 5SOP_rf3
show sticks, select_5sop
# 5SOP is rendered as loop (no secondary-structure cartoon)
cartoon loop

set_view (\
     0.410563320,    0.787584782,   -0.459498316,\
    -0.818619251,    0.096413545,   -0.566184938,\
    -0.401618332,    0.608609736,    0.684318006,\
    -0.000029873,   -0.000238253,  -47.701404572,\
    23.111997604,   21.226135254,   12.385825157,\
    38.199062347,   57.225769043,  -20.000000000 )

# GT + density mesh
enable 5SOP_isomesh
ray 1680, 1262
png ./5SOP_density.png, dpi=300
disable 5SOP_isomesh

# boltz2 X-ray
enable 5SOP_boltz2XRAY_baseline
ray 1680, 1262
png ./5SOP_boltz2XRAY_baseline.png, dpi=300
disable 5SOP_boltz2XRAY_baseline
enable 5SOP_boltz2XRAY_guided
ray 1680, 1262
png ./5SOP_boltz2XRAY_guided.png, dpi=300
disable 5SOP_boltz2XRAY_guided

# boltz2 MD
enable 5SOP_boltz2MD_baseline
ray 1680, 1262
png ./5SOP_boltz2MD_baseline.png, dpi=300
disable 5SOP_boltz2MD_baseline
enable 5SOP_boltz2MD_guided
ray 1680, 1262
png ./5SOP_boltz2MD_guided.png, dpi=300
disable 5SOP_boltz2MD_guided

# protenix
enable 5SOP_protenix_baseline
ray 1680, 1262
png ./5SOP_protenix_baseline.png, dpi=300
disable 5SOP_protenix_baseline
enable 5SOP_protenix_guided
ray 1680, 1262
png ./5SOP_protenix_guided.png, dpi=300
disable 5SOP_protenix_guided

# rf3
enable 5SOP_rf3_baseline
ray 1680, 1262
png ./5SOP_rf3_baseline.png, dpi=300
disable 5SOP_rf3_baseline
enable 5SOP_rf3_guided
ray 1680, 1262
png ./5SOP_rf3_guided.png, dpi=300
disable 5SOP_rf3_guided

# Re-enable everything for interactive inspection
enable all
