import hashlib
import tempfile
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from baseline import (
    baseline_ecdf,
    load_baseline_snapshot,
    make_baseline_interactive,
    make_baseline_static,
    prepare_baseline,
)
from loaders import RSCC_MAX_RMSD_SOURCES
from plot_distributions import make_combined_baseline_distributions
from plot_interaction import select_observations


def frames():
    common = pd.DataFrame({
        "protein": ["p1", "p1", "p2", "p3"],
        "model_name": ["protenix"] * 4,
        "selection": ["r1", "r2", "r1", "r1"],
        "selection_name": ["r1", "r2", "r1", "r1"],
        "source_file": ["input.csv"] * 4,
        "source_record": [1, 2, 3, 4],
    })
    return common.assign(rscc=[0.2, 0.8, 0.4, np.nan]), common.assign(
        clustering_score=[0, 0, 0.8, np.nan]
    )


class BaselineTest(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_population_means_missingness_and_contributors(self):
        rscc, lddt = frames()
        points, records, coverage = prepare_baseline(rscc, lddt, ["protenix"])
        rscc_points = points.loc[points.metric == "RSCC"]
        self.assertEqual(rscc_points.protein.tolist(), ["p1", "p2"])
        np.testing.assert_allclose(rscc_points.value, [0.5, 0.4])
        first = rscc_points.iloc[0]
        self.assertEqual(records.loc[records.point_id == first.point_id].source_record.tolist(), [1, 2])
        self.assertEqual(coverage.nonfinite_source_values.tolist(), [1, 1])
        regions, _, _ = prepare_baseline(rscc, lddt, ["protenix"], "region")
        self.assertEqual(len(regions.loc[regions.metric == "RSCC"]), 3)
        pd.testing.assert_frame_equal(rscc, frames()[0])

    def test_ecdf_ties_and_renderer_values_match_existing_helper(self):
        rscc, lddt = frames()
        points, _, _ = prepare_baseline(rscc, lddt, ["protenix"])
        curves = baseline_ecdf(points)
        static = make_baseline_static(curves)
        interactive = make_baseline_interactive(curves)
        reference = make_combined_baseline_distributions(
            rscc.groupby(["protein", "model_name"], as_index=False).rscc.mean(),
            lddt, cumulative=True, cumulative_style="ecdf",
        )
        for i, axis in enumerate(static.axes):
            x, y = axis.lines[0].get_data()
            np.testing.assert_array_equal(x, interactive.data[i].x)
            np.testing.assert_array_equal(y, interactive.data[i].y)
            ref_x, ref_y = reference.axes[i].lines[0].get_data()
            np.testing.assert_array_equal(x[1:], ref_x[1:])
            np.testing.assert_array_equal(y, ref_y)
            self.assertTrue(np.isfinite(x).all())
            self.assertEqual(y[-1], 1)
        cluster = curves.loc[curves.metric == "Clustering score"]
        np.testing.assert_allclose(cluster.y, [0, 1/3, 2/3, 1])
        self.assertEqual(cluster.point_id.iloc[0], "")

    def test_no_models_and_invalid_population(self):
        rscc, lddt = frames()
        points, records, coverage = prepare_baseline(rscc, lddt, [])
        self.assertTrue(baseline_ecdf(points).empty)
        self.assertTrue(records.empty)
        self.assertTrue(coverage.empty)
        with self.assertRaises(ValueError):
            prepare_baseline(rscc, lddt, ["protenix"], "unknown")

    def test_legend_contains_variants_and_groups_both_panels(self):
        rscc, lddt = frames()
        rscc = pd.concat([rscc, rscc.assign(model_name="protenix-mini")])
        lddt = pd.concat([lddt, lddt.assign(model_name="protenix-mini")])
        points, _, _ = prepare_baseline(rscc, lddt, ["protenix", "protenix-mini"])
        fig = make_baseline_interactive(baseline_ecdf(points))
        self.assertEqual(fig.layout.legend.groupclick, "togglegroup")
        for trace in fig.data:
            self.assertEqual(
                trace.visible, "legendonly" if trace.legendgroup == "protenix-mini" else True
            )
        self.assertEqual(sum(bool(t.showlegend) for t in fig.data), 2)
        self.assertEqual([t.legendgroup for t in fig.data].count("protenix-mini"), 2)

    def test_native_selection_uses_rendered_ids_even_after_table_sorting(self):
        points, _, _ = prepare_baseline(*frames(), ["protenix"])
        figure = make_baseline_interactive(baseline_ecdf(points))
        ids = [[row[0] for row in trace.customdata] for trace in figure.data]
        reordered = points.sort_values("value", ascending=False)
        for payload in [
            [{"curveNumber": 0, "pointIndex": 1}],
            [{"curveNumber": 0, "pointNumber": 1}],
            [{"customdata": [ids[0][1]]}],
        ]:
            selected = select_observations(reordered, ids, payload)
            self.assertEqual(selected.point_id.tolist(), [ids[0][1]])
        self.assertTrue(select_observations(
            reordered, ids, [{"curveNumber": 0, "pointIndex": 0}]
        ).empty)
        pd.testing.assert_frame_equal(select_observations(reordered, ids, []), reordered)

    def test_chain_selection_is_excluded_only_from_region_population(self):
        rscc, lddt = frames()
        rscc.loc[0, "selection_name"] = "all altloc"
        lddt.loc[0, "selection_name"] = "all altloc"
        _, _, coverage = prepare_baseline(rscc, lddt, ["protenix"])
        self.assertEqual(coverage.excluded_chain_rows.tolist(), [0, 1])
        _, _, coverage = prepare_baseline(rscc, lddt, ["protenix"], "region")
        self.assertEqual(coverage.excluded_chain_rows.tolist(), [1, 1])

    def test_snapshot_filters_aliases_and_never_changes_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for path, override in RSCC_MAX_RMSD_SOURCES:
                name = Path(path).name
                model = "rf3" if "rf3" in name else "protenix"
                rows = pd.DataFrame({
                    "protein": ["p1"] * 3,
                    "model": [model] * 3,
                    "method": ["X-RAY"] * 3,
                    "ensemble_size": [8] * 3,
                    "guidance_weight": [0, 0.1, 0],
                    "altloc_occupancies": [
                        "{'A': 0.5, 'B': 0.5}", "{'A': 0.5, 'B': 0.5}",
                        "{'A': 0.25, 'B': 0.75}",
                    ],
                    "rscc": [0.7, 0.9, 0.1],
                    "selection": ["chain A and resi 1-1"] * 3,
                    "occupancies": ["[0.25 0.75]"] * 3,
                    "avg_silhouette": [0.8] * 3,
                    "avg_silhouette_to_ref": [0.8] * 3,
                })
                rows.to_csv(root / name, index=False)
                lddt_name = name.replace("rscc_results_max_rmsd", "lddt_results").replace(
                    "rscc_results", "lddt_results"
                )
                rows.to_csv(root / lddt_name, index=False)
            before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir()}
            rscc, lddt, manifest = load_baseline_snapshot(root)
            self.assertEqual(len(rscc), 5)
            self.assertEqual(len(lddt), 5)
            self.assertEqual(set(rscc.model_name), {
                "protenix", "protenix-mini", "protenix-tiny", "rf3", "rf3-chiral-off",
            })
            np.testing.assert_allclose(lddt.clustering_score, 0.4)
            self.assertEqual(rscc.source_record.tolist(), [1] * 5)
            self.assertEqual(len(manifest["source_files"]), 10)
            self.assertEqual(before, {
                p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir()
            })
            (root / "lddt_results_rf3_fixed.csv").unlink()
            with self.assertRaises(FileNotFoundError):
                load_baseline_snapshot(root)


if __name__ == "__main__":
    unittest.main()
