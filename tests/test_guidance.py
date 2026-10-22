import hashlib
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.scale import SymmetricalLogTransform
import numpy as np
import pandas as pd

from guidance import (
    load_guidance_snapshot, make_guidance_interactive, make_guidance_static, prepare_guidance, symlog,
)
from loaders import RSCC_MAX_RMSD_SOURCES, add_simplified_metadata_columns
from plot_interaction import select_observations


def frames():
    rows = [
        {"protein": protein, "model_name": "protenix", "guidance_weight": gw,
         "rscc": 0.5 + (i + 1) * gw,
         "clashscore": 10 + 10 * (i + 1) * gw,
         "outlier_fraction_bond_lengths": 0.01 + (i + 1) * gw / 10,
         "outlier_fraction_bond_angles": 0.02 + (i + 1) * gw / 10}
        for i, protein in enumerate(["p1", "p2", "p3", "p4"])
        for gw in [0, 0.1, 0.2]
    ]
    df = add_simplified_metadata_columns(pd.DataFrame(rows))
    return {"rscc": df.copy(), "clashes": df.copy(), "bonds": df.copy()}


class GuidanceTest(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_deltas_then_medians_and_iqr(self):
        proteins, summary, coverage = prepare_guidance(frames())
        self.assertEqual(len(proteins), 8)
        self.assertTrue(summary.plotted.all())
        row = summary.iloc[0]
        np.testing.assert_allclose([row.x_q1, row.x_med, row.x_q3], [0.175, 0.25, 0.325])
        np.testing.assert_allclose([row.y_q1, row.y_med, row.y_q3], [1.75, 2.5, 3.25])
        self.assertEqual(row.n, 4)
        self.assertTrue((coverage.missing_baseline == 0).all())
        self.assertEqual(summary.marker_area.tolist(), [30, 220])
        self.assertEqual(summary.chosen_weight.tolist(), [True, False])

    def test_missing_baseline_and_marginal_missingness_are_visible(self):
        data = frames()
        data["rscc"] = data["rscc"].loc[
            ~((data["rscc"].protein == "p4") & (data["rscc"].guidance_weight == 0))
        ]
        data["bonds"].loc[
            (data["bonds"].protein == "p1") & (data["bonds"].guidance_weight == 0.1),
            "outlier_fraction_bond_lengths",
        ] = np.nan
        proteins, summary, coverage = prepare_guidance(data)
        self.assertEqual(summary.n.tolist(), [4, 4])
        self.assertEqual(summary.n_rscc.tolist(), [3, 3])
        self.assertEqual(summary.incomplete_geometry.tolist(), [1, 0])
        self.assertEqual(coverage.loc[coverage.metric == "rscc"].missing_baseline.tolist(), [1, 1])
        row = proteins.loc[(proteins.protein == "p1") & (proteins.guidance_weight == 0.1)].iloc[0]
        self.assertAlmostEqual(row.geometry_penalty, 1)
        self.assertEqual(row.geometry_terms, 2)

    def test_low_coverage_and_unknown_guidance(self):
        data = {key: value.loc[value.protein.isin(["p1", "p2"])] for key, value in frames().items()}
        _, summary, _ = prepare_guidance(data)
        self.assertFalse(summary.plotted.any())
        data["rscc"] = data["rscc"].copy()
        data["rscc"].loc[data["rscc"].index[0], "guidance_strength"] = None
        with self.assertRaises(ValueError):
            prepare_guidance(data)

    def test_symlog_matches_reference_at_zero_negative_and_positive_values(self):
        values = np.array([-100, -1, -0.1, 0, 0.1, 1, 2, 10, 1000])
        expected = SymmetricalLogTransform(base=10, linthresh=1, linscale=1).transform(values)
        np.testing.assert_allclose(symlog(values), expected)

    def test_both_renderers_use_final_coordinates_and_asymmetric_errors(self):
        _, summary, _ = prepare_guidance(frames())
        static = make_guidance_static(summary)
        interactive = make_guidance_interactive(summary)
        trace = interactive.data[0]
        np.testing.assert_array_equal(static.axes[0].lines[0].get_xdata(), trace.x)
        np.testing.assert_array_equal(static.axes[0].lines[0].get_ydata(), trace.y)
        np.testing.assert_allclose(trace.error_y.array[1:], summary.y_q3_plot - summary.y_med_plot)
        np.testing.assert_allclose(trace.error_y.arrayminus[1:], summary.y_med_plot - summary.y_q1_plot)
        np.testing.assert_allclose([row[2] for row in trace.customdata[1:]], summary.y_med)
        ids = [[row[0] for row in trace.customdata]]
        picked = select_observations(summary, ids, [{"curveNumber": 0, "pointIndex": 1}])
        self.assertEqual(picked.point_id.tolist(), [summary.iloc[0].point_id])

    def test_legend_keeps_variants_available_but_initially_hidden(self):
        _, summary, _ = prepare_guidance(frames())
        variant = summary.assign(model_name="protenix-mini")
        figure = make_guidance_interactive(pd.concat([summary, variant]))
        self.assertEqual([t.visible for t in figure.data], [True, "legendonly"])

    def test_snapshot_overrides_member_join_and_source_immutability(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for source, override in RSCC_MAX_RMSD_SOURCES:
                name = Path(source).name
                suffix = name.removeprefix("rscc_results_max_rmsd").removeprefix("rscc_results")
                model = "rf3" if "rf3" in name else "protenix"
                pd.DataFrame({
                    "protein": ["1abc", "1abc"], "model": [model, model], "method": ["", ""],
                    "ensemble_size": [8, 8], "guidance_weight": [0.0, 0.1],
                    "altloc_occupancies": ["{'A': 0.5, 'B': 0.5}"] * 2,
                    "rscc": [0.5, 0.6], "selection": ["chain A and resi 1-1"] * 2,
                }).to_csv(root / name, index=False)
                directory = model + "_X-RAY_DIFFRACTION" if override else model
                members = pd.DataFrame([
                    {"model": f"/data/1ABC_0.5occA_0.5occB/{directory}/pure_guidance/"
                              f"ens8_gw{gw}/refined-patched.cif",
                     "model_n": member, "outlier_fraction": 0.01 + gw}
                    for gw in [0.0, 0.1] for member in [0, 1]
                ])
                members.to_csv(root / ("bond_angle_violation_fractions" + suffix), index=False)
                members.iloc[:-1].to_csv(
                    root / ("bond_length_violation_fractions" + suffix), index=False
                )
                members.assign(
                    model_name=members.model.str.replace("refined-patched.cif", "nonan.cif", regex=False),
                    model_id=members.model_n + 1, clashscore=10.0,
                )[["model_name", "model_id", "clashscore"]].to_csv(
                    root / ("clashscore_metrics" + suffix), index=False
                )
            before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir()}
            data, manifest = load_guidance_snapshot(root)
            self.assertEqual(len(manifest["source_files"]), 20)
            self.assertEqual(len(data["bonds"]), 15)
            self.assertTrue(all(j["angle_only_members"] == 1 for j in manifest["bond_member_joins_all_occupancies"]))
            self.assertEqual(set(data["clashes"].model_name), {
                "protenix", "rf3", "protenix-mini", "protenix-tiny", "rf3-chiral-off",
            })
            self.assertEqual(before, {
                p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir()
            })
            path = root / "clashscore_metrics_rf3_fixed.csv"
            clash = pd.read_csv(path)
            pd.concat([clash, clash.iloc[:1]]).to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "duplicate clashscore"):
                load_guidance_snapshot(root)


if __name__ == "__main__":
    unittest.main()
