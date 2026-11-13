import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from fixedpoint_fir_lab.experiment import canonical_hash, make_fixtures, metrics, run_experiment, save_experiment


class ExperimentTests(unittest.TestCase):
    def test_deterministic_fixtures_and_invalid_parameters(self):
        self.assertEqual(canonical_hash(make_fixtures()), canonical_hash(make_fixtures()))
        self.assertNotEqual(canonical_hash(make_fixtures(seed=8)), canonical_hash(make_fixtures(seed=9)))
        for args in ({"seed": -1}, {"seed": True}, {"samples": 31}):
            with self.subTest(args=args), self.assertRaises(ValueError):
                make_fixtures(**args)

    def test_metrics_distinguish_rms_and_undefined_normalization(self):
        result = metrics([1, 1], [0, 0])
        self.assertEqual(result["rmse_uv"], 1)
        self.assertIsNone(result["nrmse_percent_of_reference_rms"])
        self.assertEqual(metrics([2, 2], [1, 1])["nrmse_percent_of_reference_rms"], 100)
        with self.assertRaises(ValueError):
            metrics([1], [1, 2])

    def test_quantization_stage_comparisons_are_in_uv(self):
        report, fixture = run_experiment(samples=64)
        self.assertEqual(report["fixtures"]["canonical_sha256"], canonical_hash(fixture))
        self.assertTrue(report["accumulator"]["all_direct_and_paired_partial_sums_safe_in_int32"])
        for record in report["cases"]:
            self.assertEqual(record["output_samples"], 94)
            for variant, scale in zip(record["variants"], (4, 16, 64)):
                self.assertEqual(variant["output_units_per_uv"], scale)
                self.assertTrue(variant["direct_symmetric_bitexact"])
                self.assertLess(variant["output_quantization_incremental"]["max_abs_error_uv"], 1 / scale)
        small = report["cases"][1]
        self.assertEqual(small["input_quantization_incremental"]["rmse_uv"], 0)

    def test_cli_writes_parseable_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "demo.json"
            completed = subprocess.run(
                [sys.executable, "-B", "-m", "fixedpoint_fir_lab", "experiment", "--output", str(target), "--samples", "64"],
                check=True, capture_output=True, text=True, timeout=20,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            self.assertTrue(json.loads(completed.stdout)["complete"])
            report = json.loads(target.read_text(encoding="utf-8"))
            fixture = json.loads((target.parent / "synthetic-fixtures.json").read_text(encoding="utf-8"))
            self.assertEqual(report["fixtures"]["canonical_sha256"], canonical_hash(fixture))
            self.assertFalse((target.parent / "precision.svg").exists())
            self.assertFalse(any(Path(name).is_absolute() for name in report["sources_sha256"]))
            self.assertNotIn(directory, target.read_text(encoding="utf-8"))

    def test_artifact_filename_collision_is_rejected(self):
        with self.assertRaises(ValueError):
            save_experiment("synthetic-fixtures.json")
        with self.assertRaises(ValueError):
            save_experiment("precision.svg")


if __name__ == "__main__":
    unittest.main()
