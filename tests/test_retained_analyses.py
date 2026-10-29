import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from loaders import get_rscc_nan_proteins, merge_classification_by_containment
from plot_scatter import make_min_rmsd_scatter
from plot_training_bias import make_training_bias_no_cluster_strip
from transforms import build_training_bias_df
from transforms import filter_by_model_group


def notebook(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).parents[1] / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RetainedAnalysesTest(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_classification_keeps_narrowest_match_and_drops_unmatched_and_unparseable(self):
        rscc = pd.DataFrame({
            "protein": ["p", "p", "p", "p", "domain"],
            "selection": ["chain A and resi 1-3", "chain A and resi 2-2",
                          "chain A and resi 8-10", "all altloc", "chain B and resi 100-103"],
            "rscc": [.5, .6, .7, .8, .9],
        })
        classes = pd.DataFrame({
            "protein": ["p", "p", "domain"], "chain": ["A", "A", "B"],
            "selection": ["wide", "narrow", "domain"],
            "start_res": [1, 1, 100], "end_res": [5, 3, 105],
            "classification": ["large_loop", "small_loop", "domain_shift"],
        })
        result = merge_classification_by_containment(rscc, classes)
        self.assertEqual(dict(zip(result.selection, result.classification)), {
            "chain A and resi 1-3": "small_loop",
            "chain A and resi 2-2": "small_loop",
            "chain B and resi 100-103": "domain_shift",
        })
        self.assertEqual(set(result.rscc), {.5, .6, .9})

    def test_nan_blacklist_includes_rf3_but_not_its_alternative_variant(self):
        rscc = pd.DataFrame({
            "protein": ["rf3_only", "boltz_only", "alternative_only", "finite"],
            "model_name": ["rf3", "boltz2MD", "rf3-chiral-off", "protenix"],
            "rscc": [np.nan, np.nan, np.nan, .5],
        })
        with patch("loaders.load_rscc_max_rmsd", return_value=rscc):
            self.assertEqual(get_rscc_nan_proteins(), {"rf3_only", "boltz_only"})
            self.assertEqual(get_rscc_nan_proteins(("boltz2MD",)), {"boltz_only"})

    def test_training_bias_uses_second_occupancy_not_minimum_and_keeps_cluster_states(self):
        lddt = pd.DataFrame({
            "protein": ["empty", "unassigned", "one_sided", "mixed"], "selection": ["r"] * 4,
            "model_name": ["protenix"] * 4, "guidance_weight": [.1] * 4,
            "occupancies": ["[0.25 0.75]", "[0.75 0.25]", "[0.125 0.875]", "[0.5 0.5]"],
        })
        clusters = pd.DataFrame({
            "protein": ["empty", "unassigned", "one_sided", "mixed"], "selection": ["r"] * 4,
            "n_assigned_to_A": [0, 0, 3, 3], "n_assigned_to_B": [0, 0, 0, 2],
            "n_cluster_evaluated": [0, 10, 3, 5],
            "frac_a_cluster": [np.nan, np.nan, 1., .6], "x_one_minus_fa": [np.nan, np.nan, 0., .4],
        })
        cutoffs = pd.DataFrame({
            "protein": ["empty", "unassigned", "one_sided", "mixed"],
            "cutoff_status": ["Pre", "Pre", "Post", "Pre"],
        })
        result = build_training_bias_df(
            lddt, clusters, cutoffs, {"protenix": .1}, verbose=False,
        ).set_index("protein")
        np.testing.assert_allclose(result.b_occupancy, [.75, .25, .875, .5])
        self.assertTrue(np.isnan(result.loc["empty", "x_one_minus_fa"]))
        self.assertEqual(result.loc["one_sided", "n_assigned_to_B"], 0)
        self.assertEqual(result.loc["one_sided", "cutoff_status"], "Post")
        self.assertEqual(result.loc["mixed", "x_one_minus_fa"], .4)
        # Preserve Figure 32's current A=B=0 population until scientific adjudication.
        figure = make_training_bias_no_cluster_strip(
            result.reset_index(), result.reset_index(), model_order=["protenix", "rf3"],
        )
        ys = pd.DataFrame(figure.layout.meta["observations"]).y
        np.testing.assert_allclose(sorted(ys), [.25, .25, .75, .75])

    def test_per_model_rmsd_still_draws_unaggregated_values_in_both_panels(self):
        points = pd.DataFrame({
            "min_rmsd_to_A": [1., 2.], "min_rmsd_to_B": [2., 1.], "label": ["r1", "r2"],
            "protein": ["p1", "p2"], "selection": ["r1", "r2"],
        })
        figure = make_min_rmsd_scatter(points, title="Per-model", pinned_highlights=())
        copies = [t for t in figure.data if t.meta and t.meta.get("observations")]
        self.assertEqual(len(copies), 2)
        for trace in copies:
            np.testing.assert_array_equal(np.column_stack([trace.x, trace.y]), [[1., 2.], [2., 1.]])

    def test_training_inputs_apply_blacklist_to_every_model_group_and_metadata(self):
        cell = notebook("07_training_tradeoffs").load_training_inputs
        metadata = pd.DataFrame({"protein": ["drop", "keep"]})

        def load_lddt(path, model_name_override=None):
            model = model_name_override or ("rf3" if "rf3_fixed" in path else "protenix")
            return metadata.assign(model_name=model)

        for group, expected_models in [
            ("main", {"protenix", "rf3"}),
            ("protenix", {"protenix", "protenix-mini", "protenix-tiny"}),
            ("rf3", {"rf3", "rf3-chiral-off"}),
        ]:
            for excluded in (set(), {"drop"}):
                _, definitions = cell.run(
                    filter_by_model_group=filter_by_model_group,
                    load_deposition_cutoffs=lambda: metadata.copy(),
                    load_training_cluster_counts=lambda: metadata.copy(),
                    load_lddt=load_lddt, pd=pd, model_group=SimpleNamespace(value=group),
                    nan_rscc_proteins=excluded,
                )
                for key in ("cluster_df", "cutoffs_df", "lddt_for_bias"):
                    self.assertEqual(set(definitions[key].protein), {"drop", "keep"} - excluded)
                self.assertEqual(set(definitions["lddt_for_bias"].model_name), expected_models)

    def test_budget_grid_wins_over_default_without_rowwise_backfill(self):
        def row(model, selection, rscc, recycles=10):
            return {"protein": "5mhx", "model_name": model, "guidance_weight": .1,
                    "n_recycles": recycles, "n_diffusion_steps": 200,
                    "selection": selection, "rscc": rscc}

        grid = pd.DataFrame([row("protenix", "r1", .9), row("protenix", "r1", .7, recycles=1)])
        default = pd.DataFrame([
            row("protenix", "r1", .2), row("protenix", "r2", .3), row("rf3", "r1", .4),
        ])
        _, definitions = notebook("08_sampling_budgets").combine_budget_grid.run(
            default_df=default, grid_df=grid, pd=pd,
        )
        combined = definitions["combined_df"]
        self.assertEqual(len(combined), 3)
        self.assertEqual(set(combined.rscc), {.9, .7, .4})
        self.assertNotIn("r2", combined.selection.tolist())


if __name__ == "__main__":
    unittest.main()
