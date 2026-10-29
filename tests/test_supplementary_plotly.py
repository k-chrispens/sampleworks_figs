from io import BytesIO
import json
import unittest
from unittest.mock import patch
from zipfile import ZipFile

import marimo as mo
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from guidance import symlog
from plot_budgets import make_budget_figure
from plot_classified import make_by_class_combined_grid
from plot_interaction import highlight_observations, highlight_table_rows, select_observations
from plot_scatter import make_rscc_scatter, make_clustering_scatter, rscc_residual_label_override
from plot_tradeoff import make_spiral_figure
from plot_training_bias import (
    _training_bias_binned_means, make_training_bias_scatter,
    make_training_bias_calibration_overlay, make_training_bias_small_multiples_v2,
    make_training_bias_one_sided_strip, make_training_bias_no_cluster_strip,
    make_training_bias_calibration_explainer, make_training_bias_calibration_overlay_legend,
)
from supplementary_plot import linked_plot, observations, reconstruction_bundle, inspection_tables, outline_rows


def points(figure):
    return pd.DataFrame(figure.layout.meta["observations"])


def summaries(figure, kind):
    frame = pd.DataFrame(figure.layout.meta["summaries"])
    return frame.loc[frame.kind == kind]


def training():
    return pd.DataFrame({
        "protein": [f"p{i}" for i in range(8)], "selection": ["r"] * 8,
        "model_name": ["protenix"] * 8, "guidance_weight": [.1] * 8,
        "n_assigned_to_A": [0, 1, 2, 3, 0, 0, 2, 0],
        "n_assigned_to_B": [4, 3, 2, 1, 0, 0, 0, 0],
        "n_cluster_evaluated": [5, 5, 6, 6, 0, 10, 2, 0],
        "b_occupancy": [.8, .4, .6, .2, .25, .75, .1, .9],
        "cutoff_status": ["Pre", "Post", "Pre", "Post", "Pre", "Post", "Pre", "Pre"],
    })


class SupplementaryPlotlyTest(unittest.TestCase):
    def test_residual_values_densities_and_mirrored_ring_sets(self):
        data = pd.DataFrame({
            "protein": [f"p{i}" for i in range(8)], "selection": ["r"] * 8,
            "model_name": ["protenix"] * 8, "label": [f"P{i} r" for i in range(8)],
            "rscc_baseline": [.1, .2, .3, .4, .5, .6, .7, .8],
            "rscc_optimal": [.9, .8, .6, .2, .3, .7, .5, .8],
        })
        figure = make_rscc_scatter(data, "protenix", .1)
        np.testing.assert_allclose(points(figure).y, data.rscc_optimal - data.rscc_baseline)
        kde = gaussian_kde(data.rscc_baseline)
        kde.set_bandwidth(kde.factor * .5)
        bandwidth = np.sqrt(kde.covariance.squeeze())
        support = np.linspace(.1 - 3 * bandwidth, .8 + 3 * bandwidth, 200)
        density = summaries(figure, "density").query("metric == 'rscc_baseline'")
        np.testing.assert_allclose(density.y, support, atol=1e-12)
        np.testing.assert_allclose(density.x, kde(support), atol=1e-12)
        clustering = data.rename(columns={"rscc_baseline": "clustering_score_baseline",
                                         "rscc_optimal": "clustering_score_optimal"})
        clustering["clustering_score_optimal"] = 1 - clustering.clustering_score_optimal
        mirrored = make_clustering_scatter(clustering, "protenix", .1,
                                           residual_label_override=rscc_residual_label_override(data))
        self.assertEqual(points(figure).highlight_colors.tolist(), points(mirrored).highlight_colors.tolist())
        self.assertFalse(any(t.mode and "text" in t.mode for t in figure.data))

    def test_class_grid_shadows_and_guides_are_not_observations(self):
        data = pd.DataFrame({
            "protein": ["p", "q"], "selection": ["r", "r"], "model_name": ["protenix"] * 2,
            "classification": ["small_loop", "domain_shift"],
            "rscc_baseline": [.2, .8], "rscc_optimal": [.6, .7],
        })
        figure = make_by_class_combined_grid(data, model_order=["protenix"])
        foreground = [t for t in figure.data if t.meta and t.meta.get("observations")]
        self.assertEqual(sum(len(t.x) for t in foreground), 2)
        np.testing.assert_allclose(points(figure).y, [.4, -.1])
        ids = [[r[0] for r in t.customdata] for t in figure.data]
        self.assertTrue(select_observations(points(figure), ids, [{"curveNumber": 0, "pointIndex": 0}]).empty)

    def test_trajectories_keep_raw_values_and_transform_only_display(self):
        frame = pd.DataFrame({
            "protein": ["p"] * 3, "model_name": ["protenix"] * 3,
            "guidance_weight": [.2, 0, .1], "delta_rscc": [.4, 0, .2],
            "geometry_penalty": [100, 0, -1],
        })
        figure = make_spiral_figure(frame, "protenix", xlim=(-.1, .5), ylim=(-2, 200))
        prepared = points(figure)
        np.testing.assert_allclose(prepared.geometry_penalty, [0, -1, 100])
        np.testing.assert_allclose(prepared.y, symlog([0, -1, 100]))
        self.assertEqual(prepared.guidance_weight.tolist(), [0, .1, .2])
        self.assertEqual(figure.data[0].x[0], 0)

    def test_quantile_bins_keep_ties_and_sample_sd(self):
        x = np.array([0, 0, .5, .5, 1, 1])
        y = np.array([0, .2, .4, .6, .8, 1])
        centers, means, sd = _training_bias_binned_means(x, y, n_bins=3, min_count=2)
        np.testing.assert_allclose(centers, [0, .5, 1])
        np.testing.assert_allclose(means, [.1, .5, .9])
        np.testing.assert_allclose(sd, [np.sqrt(.02)] * 3)

    def test_calibration_summaries_ignore_jitter_and_filter_inclusively(self):
        data = training()
        figures = [make_training_bias_calibration_overlay(data, data, model_order=["protenix"],
                                                          jitter=j, n_bins=2, min_bin_count=2)
                   for j in (0, .01)]
        pd.testing.assert_frame_equal(summaries(figures[0], "bin"), summaries(figures[1], "bin"))
        pd.testing.assert_frame_equal(summaries(figures[0], "regression"),
                                      summaries(figures[1], "regression"))
        stats = summaries(figures[0], "regression")
        raw = data.loc[(data.n_assigned_to_A + data.n_assigned_to_B) > 0]
        x = raw.n_assigned_to_B / (raw.n_assigned_to_A + raw.n_assigned_to_B)
        np.testing.assert_allclose(stats.slope, np.polyfit(x, raw.b_occupancy, 1)[0])
        histogram = summaries(figures[0], "histogram")
        self.assertEqual(histogram.x.sum(), len(raw) * 2)
        low = make_training_bias_calibration_overlay(data, data, model_order=["protenix"], max_cluster_size=5)
        self.assertEqual(set(points(low).protein), {"p0", "p1", "p6"})
        post = make_training_bias_calibration_overlay(data, data, model_order=["protenix"], cutoff_status="Post")
        self.assertEqual(set(points(post).protein), {"p1", "p3"})

    def test_every_training_layout_handles_observations_and_empty_models(self):
        data = training()
        for make in (make_training_bias_scatter, make_training_bias_calibration_overlay,
                     make_training_bias_small_multiples_v2, make_training_bias_no_cluster_strip,
                     make_training_bias_one_sided_strip):
            figure = make(data, data, model_order=["protenix", "rf3"])
            self.assertTrue(points(figure).point_id.is_unique)
            empty = make(data.iloc[:0], data.iloc[:0], model_order=["protenix"])
            self.assertTrue(points(empty).empty)
        empty_cluster = points(make_training_bias_no_cluster_strip(data, data, model_order=["protenix"]))
        self.assertEqual(set(empty_cluster.protein), {"p4", "p5", "p7"})
        sided = points(make_training_bias_one_sided_strip(data, data, model_order=["protenix"]))
        self.assertEqual(set(sided.protein), {"p0", "p6"})
        for make in (make_training_bias_calibration_explainer, make_training_bias_calibration_overlay_legend):
            self.assertTrue(points(make()).empty)
        legend = make_training_bias_calibration_overlay_legend(model_order=["protenix"])
        self.assertEqual([trace.name for trace in legend.data], ["Baseline", "Protenix"])
        self.assertTrue(all(trace.showlegend and list(trace.x) == [None] for trace in legend.data))

    def test_budgets_sample_sd_held_axes_and_contributors(self):
        data = pd.DataFrame([
            dict(model_name="protenix", protein="p", selection=f"r{i}", rscc=value,
                 n_recycles=rec, n_diffusion_steps=steps)
            for rec, steps, values in [(1, 200, [.2, .6]), (10, 200, [.5, .9]), (10, 50, [.1, .3])]
            for i, value in enumerate(values)
        ])
        figure = make_budget_figure(data, ["protenix"], default_recycles={"protenix": 10})
        prepared = points(figure)
        self.assertEqual(prepared.n.tolist(), [2] * 4)
        np.testing.assert_allclose(prepared["mean"], [.4, .7, .2, .7])
        np.testing.assert_allclose(prepared.sd, [np.sqrt(.08), np.sqrt(.08), np.sqrt(.02), np.sqrt(.08)])
        contributors = pd.DataFrame(figure.layout.meta["contributors"])
        self.assertEqual(len(contributors), 8)
        for point_id, group in contributors.groupby("point_id"):
            actual = prepared.set_index("point_id").loc[point_id]
            self.assertAlmostEqual(actual["mean"], group.rscc.mean())
        ensemble = data.rename(columns={"n_recycles": "ensemble_size"}).drop(columns="n_diffusion_steps")
        ensemble = ensemble.iloc[:4]
        result = points(make_budget_figure(ensemble, ["protenix"]))
        np.testing.assert_allclose(result.x, np.log2(result.budget_value))

    def test_linked_highlighting_preserves_view_and_fixed_export_data(self):
        data = training()
        base = make_training_bias_calibration_overlay(data, data, model_order=["protenix"])
        prepared = points(base)
        panel = linked_plot(base, mo, lambda _: None, "test")
        panel.widget.layout.xaxis.range = [.1, .9]
        panel.widget.data[0].visible = "legendonly"
        before = base.to_json()
        highlight_table_rows(panel.widget, base, prepared.iloc[[0]])
        self.assertEqual(base.to_json(), before)
        self.assertEqual(panel.widget.layout.xaxis.range, (.1, .9))
        self.assertEqual(panel.widget.data[0].visible, "legendonly")
        selected = highlight_observations(base, prepared.point_id.iloc[:1])
        self.assertTrue(any(t.marker.line.width is not None for t in selected.data if t.meta))
        with patch.object(type(base), "to_image", return_value=b"<svg/>"):
            content = reconstruction_bundle(base)
        with ZipFile(BytesIO(content)) as bundle:
            exported = pd.read_csv(bundle.open("points.csv"))
            self.assertEqual(len(exported), len(prepared))
            recipe = json.loads(bundle.read("recipe.json"))
            self.assertIn("no transient", recipe["scope"])

    def test_duplicate_scientific_keys_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate plotted"):
            observations(pd.DataFrame({"protein": ["p", "p"], "selection": ["r", "r"]}))

    def test_clear_is_an_event_even_without_a_previous_lasso_and_old_events_do_not_leak(self):
        base = make_training_bias_calibration_overlay(training(), training(), model_order=["protenix"])
        panel = linked_plot(base, mo, lambda _: None, "training")
        inspection_tables({"training": panel}, {}, mo)
        original_table = panel.table
        outline_rows(panel, panel.points.iloc[[0]])
        self.assertTrue(panel.has_outlines)
        inspection_tables({"training": panel}, {
            "training": {"generation": panel.generation, "revision": "clear-1", "points": []},
        }, mo)
        self.assertFalse(panel.has_outlines)
        self.assertIsNot(panel.table, original_table)
        rebuilt = linked_plot(base, mo, lambda _: None, "training")
        with patch.object(mo.ui, "table", wraps=mo.ui.table) as table:
            inspection_tables({"training": rebuilt}, {
                "training": {"generation": panel.generation, "revision": "old-selection",
                             "points": [{"customdata": [panel.points.point_id.iloc[0]]}]},
            }, mo)
        self.assertEqual(len(table.call_args.args[0]), len(rebuilt.points))
        self.assertFalse(rebuilt.has_outlines)

    def test_native_legend_groups_include_summaries_and_uncertainty(self):
        figure = make_training_bias_calibration_overlay(training(), training(), model_order=["protenix"],
                                                        n_bins=2, min_bin_count=2)
        for group in ("protenix:Baseline", "protenix:Optimal guidance"):
            traces = [trace for trace in figure.data if trace.legendgroup == group]
            self.assertEqual(len(traces), 5)  # points, lower/upper SD, means, histogram
            self.assertEqual(sum(bool(trace.showlegend) for trace in traces), 1)


if __name__ == "__main__":
    unittest.main()
