import hashlib
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from loaders import load_min_altloc_rmsd
from min_rmsd import (
    highlighted_min_rmsd, load_min_rmsd_snapshot, make_min_rmsd_interactive,
    make_min_rmsd_static, prepare_min_rmsd,
)
from plot_interaction import select_observations
from style import OPTIMAL_GW


def frame(protein="6b8x", selection="chain A and resi 181-183"):
    return pd.DataFrame([
        {"protein": protein, "selection": selection, "selection_name": selection,
         "model_name": model, "guidance_weight": OPTIMAL_GW[model],
         "min_rmsd_to_A": a, "min_rmsd_to_B": b,
         "source_file": "input.csv", "source_record": i + 1}
        for i, (model, a, b) in enumerate([
            ("boltz2MD", 1., 2.), ("boltz2X-RAY", 2., 4.),
            ("protenix", 9., 6.), ("rf3", 12., 8.),
        ])
    ])


class MinRmsdTest(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_mean_and_median_are_across_models_not_another_minimum(self):
        data = frame()
        median, _, _ = prepare_min_rmsd(data)
        mean, _, _ = prepare_min_rmsd(data, "mean")
        self.assertEqual(median.iloc[0].min_rmsd_to_A, 5.5)
        self.assertEqual(mean.iloc[0].min_rmsd_to_A, 6.)
        self.assertEqual(mean.iloc[0].min_rmsd_to_B, 5.)
        self.assertEqual(mean.iloc[0].n_models, 4)
        pd.testing.assert_frame_equal(data, frame())

    def test_failures_excluded_as_pairs_and_missing_conditions_reported(self):
        data = frame()
        data.loc[0, "min_rmsd_to_A"] = np.nan
        points, contributions, coverage = prepare_min_rmsd(data)
        self.assertEqual(points.iloc[0].n_models, 3)
        self.assertEqual(points.iloc[0].min_rmsd_to_B, 6.)
        self.assertFalse(contributions.iloc[0].included)
        self.assertEqual(coverage.nonfinite_pairs.sum(), 1)
        data.loc[1, "guidance_weight"] = 0
        points, _, coverage = prepare_min_rmsd(data)
        self.assertFalse(points.iloc[0].plotted)
        self.assertEqual(coverage.missing_chosen_weight.sum(), 1)

    def test_invalid_identity_duplicate_and_negative_values_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            prepare_min_rmsd(pd.concat([frame(), frame().iloc[:1]]))
        for column, value in [("selection", None), ("min_rmsd_to_A", -1.)]:
            data = frame()
            data.loc[0, column] = value
            with self.assertRaises(ValueError):
                prepare_min_rmsd(data)
        with self.assertRaisesRegex(ValueError, "aggregation"):
            prepare_min_rmsd(frame(), "min")

    def test_highlight_table_and_rendered_selection_match_both_panels(self):
        blue = frame()
        blue[["min_rmsd_to_A", "min_rmsd_to_B"]] = [1., 2.]
        red = frame("7z0e", "chain A and resi 11-13")
        red[["min_rmsd_to_A", "min_rmsd_to_B"]] = [2., 1.]
        points, _, _ = prepare_min_rmsd(pd.concat([blue, red]))
        table = highlighted_min_rmsd(points)
        self.assertEqual(len(table), 2)
        self.assertTrue(table.point_id.is_unique)
        self.assertIn("zoom only", table.loc[table.protein == "7z0e", "rings"].iloc[0])
        static = make_min_rmsd_static(points)
        interactive = make_min_rmsd_interactive(points)
        np.testing.assert_allclose(static.axes[0].collections[0].get_offsets(),
                                   np.array([interactive.data[0].x, interactive.data[0].y]).T)
        for axis in static.axes:
            self.assertEqual({t.get_text() for t in axis.texts}, {"closer to A", "closer to B"})
        ids = [[r[0] for r in trace.customdata] for trace in interactive.data]
        base_traces = [i for i, t in enumerate(interactive.data) if t.name == "Regions"]
        picked = select_observations(points, ids, [
            {"curveNumber": i, "pointIndex": 0} for i in base_traces
        ])
        self.assertEqual(len(picked), 1)
        red_traces = [t for t in interactive.data if "zoom only" in t.name]
        self.assertEqual(len(red_traces), 1)
        self.assertEqual(red_traces[0].xaxis, "x2")
        self.assertEqual(interactive.layout.yaxis2.constrain, "domain")
        self.assertEqual(tuple(interactive.layout.yaxis2.range), (0, 3))

    def test_all_failed_is_an_empty_analysis_not_an_invented_value(self):
        data = frame()
        data["min_rmsd_to_A"] = np.nan
        points, contributions, coverage = prepare_min_rmsd(data)
        self.assertTrue(points.empty)
        self.assertEqual(len(contributions), 4)
        self.assertEqual(coverage.nonfinite_pairs.sum(), 4)
        self.assertTrue(highlighted_min_rmsd(points).empty)

    def test_snapshot_keeps_failure_evidence_without_changing_legacy_loader(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for suffix, model in [("boltz_protenix", "protenix"), ("rf3_fixed", "rf3")]:
                pd.DataFrame({
                    "protein": ["1abc"] * 2, "model": [model] * 2, "method": [""] * 2,
                    "ensemble_size": [8] * 2, "guidance_weight": [OPTIMAL_GW[model]] * 2,
                    "altloc_occupancies": ["{'A': 0.5, 'B': 0.5}"] * 2,
                    "selection": ["chain A and resi 1-3", "chain A and resi 4-6"],
                    "min_rmsd_to_A": [.2, np.nan], "min_rmsd_to_B": [.3, .4],
                }).to_csv(root / f"min_altloc_rmsd_results_{suffix}.csv", index=False)
            before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir()}
            data, manifest = load_min_rmsd_snapshot(root)
            self.assertEqual(len(data), 4)
            self.assertEqual(data.source_record.tolist(), [1, 2, 1, 2])
            self.assertEqual(len(manifest["source_files"]), 2)
            self.assertEqual(len(load_min_altloc_rmsd(str(root / "min_altloc_rmsd_results_rf3_fixed.csv"))), 1)
            self.assertEqual(before, {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir()})


if __name__ == "__main__":
    unittest.main()
