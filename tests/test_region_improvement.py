import hashlib
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plot_interaction import select_observations
from region_improvement import (
    HIGHLIGHTS, highlighted_regions,
    load_improvement_snapshot, make_improvement_interactive,
    make_improvement_static, prepare_density, prepare_improvement,
)
from style import OPTIMAL_GW

MODELS = ["boltz2MD", "boltz2X-RAY", "protenix", "rf3"]


def frame(baselines=(0.0, 0.4, 0.6, 0.7), guided=(0.2, 0.9, 0.65, 0.8),
          protein="1abc", selection="chain A and resi 1-3"):
    return pd.DataFrame([
        {"protein": protein, "model_name": model, "selection": selection,
         "selection_name": selection, "guidance_weight": gw, "rscc": value,
         "source_file": "input.csv", "source_record": 2 * i + offset + 1}
        for i, model in enumerate(MODELS)
        for offset, (gw, value) in enumerate([(0, baselines[i]), (OPTIMAL_GW[model], guided[i])])
    ])


class RsccImprovementTest(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_difference_of_medians_not_median_of_differences(self):
        original = frame()
        points, pairs, _ = prepare_improvement(original)
        row = points.iloc[0]
        self.assertAlmostEqual(row.baseline_value, 0.5)
        self.assertAlmostEqual(row.guided_value, 0.725)
        self.assertAlmostEqual(row.improvement, 0.225)
        self.assertAlmostEqual(pairs.delta_value.median(), 0.15)
        self.assertEqual(row.n_model_rows, 4)
        pd.testing.assert_frame_equal(original, frame())

    def test_missing_pairs_and_legacy_row_threshold_are_reported(self):
        data = frame()
        data = data.loc[~((data.model_name == "rf3") & (data.guidance_weight == 0.01))].copy()
        data.loc[(data.model_name == "protenix") & (data.guidance_weight > 0), "rscc"] = np.nan
        points, _, coverage = prepare_improvement(data)
        self.assertTrue(points.iloc[0].plotted)
        self.assertEqual(points.iloc[0].n_model_rows, 3)
        self.assertEqual(points.iloc[0].n_finite_pairs, 2)
        self.assertEqual(coverage.loc[coverage.model_name == "rf3"].missing_guided.iloc[0], 1)
        points, _, _ = prepare_improvement(data.loc[data.model_name != "protenix"])
        self.assertFalse(points.iloc[0].meets_minimum)

    def test_duplicate_pairs_and_missing_identity_are_rejected(self):
        data = frame()
        with self.assertRaises(ValueError):
            prepare_improvement(pd.concat([data, data.iloc[:1]]))
        data.loc[0, "selection"] = None
        with self.assertRaisesRegex(ValueError, "identity"):
            prepare_improvement(data)

    def test_highlight_precedence_and_no_far_candidate_backfill(self):
        rows = [
            frame([x] * 4, [y] * 4, protein=p, selection=s)
            for x, y, p, s in [
                (.1, .95, "p1", "r1"), (.2, .9, "p2", "r1"),
                (.3, .85, "6b8x", "chain A and resi 181-183"),
                (.9, .8, "p4", "r1"), (.5, .3, "p5", "r1"),
                (.1, .2, "p6", "r1"), (.8, .8, "p7", "r1"),
            ]
        ]
        points, _, _ = prepare_improvement(pd.concat(rows))
        self.assertEqual(set(points.loc[points.top].protein), {"p1", "p2", "6b8x"})
        self.assertEqual(set(points.loc[points.decreased].protein), {"p4", "p5"})
        self.assertEqual(set(points.loc[points.far].protein), {"p6", "p7"})
        curated = points.loc[points.protein == "6b8x"].iloc[0]
        self.assertTrue(curated.top and curated.curated)
        table = highlighted_regions(points)
        self.assertEqual(len(table), 7)
        self.assertTrue(table.point_id.is_unique)
        rings = table.loc[table.protein == "6b8x", "rings"].iloc[0]
        self.assertIn("Most improved (#2E7D32)", rings)
        self.assertIn("Biological examples (#1565C0)", rings)
        points.loc[points.protein == "6b8x", "plotted"] = False
        self.assertNotIn("6b8x", highlighted_regions(points).protein.tolist())

    def test_density_and_renderers_share_final_values(self):
        points, _, _ = prepare_improvement(pd.concat([
            frame(protein="p1"), frame([.2, .3, .4, .5], [.4, .5, .6, .8], protein="p2"),
            frame([.5, .6, .7, .8], [.8, .7, .8, .9], protein="p3"),
        ]))
        density, audit = prepare_density(points)
        for series, column in [("Baseline", "baseline_value"), ("Guided", "guided_value")]:
            _, axis = plt.subplots()
            sns.kdeplot(points[column], ax=axis, bw_adjust=0.5)
            curve = density.loc[density.series == series]
            np.testing.assert_allclose(curve.x, axis.lines[0].get_xdata(), rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(curve.density, axis.lines[0].get_ydata(), rtol=1e-12, atol=1e-12)
        self.assertTrue(audit.status.eq("drawn").all())
        static = make_improvement_static(points, density)
        interactive = make_improvement_interactive(points, density)
        self.assertEqual(len(static.axes[0].texts), 0)
        for collection, (_, flag, _) in zip(static.axes[0].collections[1:], HIGHLIGHTS):
            self.assertEqual(len(collection.get_offsets()), int(points[flag].sum()))
        self.assertTrue(all(trace.mode != "markers+text" for trace in interactive.data))
        self.assertEqual(interactive.data[0].marker.line.width, 0)
        np.testing.assert_allclose(static.axes[0].collections[0].get_offsets()[:, 0], interactive.data[0].x)
        np.testing.assert_allclose(static.axes[0].collections[0].get_offsets()[:, 1], interactive.data[0].y)
        ids = [[row[0] for row in trace.customdata] for trace in interactive.data]
        picked = select_observations(points, ids, [{"curveNumber": 0, "pointIndex": 0}])
        self.assertEqual(picked.point_id.tolist(), [points.iloc[0].point_id])
        self.assertTrue(select_observations(
            points, ids, [{"curveNumber": len(ids) - 1, "pointIndex": 0}]
        ).empty)

    def test_constant_density_is_explicitly_unavailable(self):
        points, _, _ = prepare_improvement(frame())
        density, audit = prepare_density(points)
        self.assertTrue(density.empty)
        self.assertTrue(audit.status.str.startswith("unavailable").all())

    def test_no_region_pairs_produces_an_explicit_empty_analysis(self):
        data = frame()
        data["selection_name"] = "all altloc"
        points, pairs, coverage = prepare_improvement(data)
        self.assertTrue(points.empty)
        self.assertTrue(pairs.empty)
        self.assertTrue(coverage.empty)
        density, audit = prepare_density(points)
        self.assertTrue(density.empty)
        self.assertTrue(audit.status.str.startswith("unavailable").all())

    def test_snapshot_uses_only_main_files_and_preserves_them(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for suffix, model in [("boltz_protenix", "protenix"), ("rf3_fixed", "rf3")]:
                pd.DataFrame({
                    "protein": ["1abc", "1abc"], "model": [model, model], "method": ["", ""],
                    "ensemble_size": [8, 8], "guidance_weight": [0.0, OPTIMAL_GW[model]],
                    "altloc_occupancies": ["{'A': 0.5, 'B': 0.5}"] * 2,
                    "rscc": [0.5, 0.7], "selection": ["chain A and resi 1-3"] * 2,
                }).to_csv(root / f"rscc_results_max_rmsd_{suffix}.csv", index=False)
            before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir()}
            data, manifest = load_improvement_snapshot(root)
            self.assertEqual(len(data), 4)
            self.assertEqual(len(manifest["source_files"]), 2)
            self.assertEqual(before, {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir()})

    def test_clustering_uses_shared_coordinates_and_retains_score_components(self):
        data = frame().rename(columns={"rscc": "clustering_score"})
        data["occupancies"] = "[0.5 0.5]"
        data["avg_silhouette_to_ref"] = data.clustering_score
        data["min_predicted_occupancy"] = .5
        data["min_map_occupancy"] = .5
        points, pairs, _ = prepare_improvement(data, "clustering_score")
        self.assertAlmostEqual(points.iloc[0].improvement, .225)
        self.assertIn("occupancies_baseline", pairs)
        density, _ = prepare_density(points)
        static = make_improvement_static(points, density, "clustering_score")
        interactive = make_improvement_interactive(points, density, "clustering_score")
        self.assertIn("Clustering score", static.axes[1].get_xlabel())
        self.assertIn("Clustering score", interactive.layout.xaxis2.title.text)
        self.assertEqual(len(static.axes[0].texts), 0)

    def test_far_candidate_count_is_three_for_clustering_and_five_for_rscc(self):
        data = pd.concat([
            frame([value] * 4, [value] * 4, protein=f"p{i}")
            for i, value in enumerate([.1, .2, .3, .4, .5, .6])
        ])
        rscc, _, _ = prepare_improvement(data)
        clustering = data.rename(columns={"rscc": "clustering_score"})
        clustering["occupancies"] = "[0.5 0.5]"
        clustering["avg_silhouette_to_ref"] = clustering.clustering_score
        clustering["min_predicted_occupancy"] = .5
        clustering["min_map_occupancy"] = .5
        scores, _, _ = prepare_improvement(clustering, "clustering_score")
        self.assertEqual(set(rscc.loc[rscc.far, "protein"]), {"p0", "p1", "p2", "p3", "p4"})
        self.assertEqual(set(scores.loc[scores.far, "protein"]), {"p0", "p1", "p2"})

    def test_clustering_snapshot_preserves_multi_entry_formula_and_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for suffix, model in [("boltz_protenix", "protenix"), ("rf3_fixed", "rf3")]:
                pd.DataFrame({
                    "protein": ["1abc"] * 3, "model": [model] * 3, "method": [""] * 3,
                    "ensemble_size": [8] * 3, "guidance_weight": [0, .1, .2],
                    "altloc_occupancies": ["{'A': 0.3, 'B': 0.7}", "{'A': 0.5, 'B': 0.5}", "{'A': 0.5, 'B': 0.5}"],
                    "occupancies": ["[0.5 0.5]", "[0.5 0.25 0.25]", ""],
                    "avg_silhouette": [.8] * 3, "avg_silhouette_to_ref": [.8] * 3,
                    "selection": ["chain A and resi 1-3"] * 3,
                }).to_csv(root / f"lddt_results_{suffix}.csv", index=False)
            before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir()}
            data, manifest = load_improvement_snapshot(root, "clustering_score")
            self.assertEqual(data.source_record.tolist(), [2, 3, 2, 3])
            np.testing.assert_allclose(data.iloc[[0, 2]].clustering_score, [.4, .4])
            self.assertTrue(data.iloc[[1, 3]].clustering_score.isna().all())
            self.assertEqual(len(manifest["source_files"]), 2)
            self.assertEqual(before, {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir()})


if __name__ == "__main__":
    unittest.main()
