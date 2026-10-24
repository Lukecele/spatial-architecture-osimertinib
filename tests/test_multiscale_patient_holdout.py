import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "extension_pilot" / "multiscale_patient_holdout_results.json"


class TestMultiscalePatientHoldout(unittest.TestCase):
    def test_pre_specified_window_is_patient_held_out(self):
        if not OUTPUT.exists():
            subprocess.run(
                [sys.executable, "src/extension_multiscale_patient_holdout.py"],
                cwd=ROOT,
                check=True,
            )
        data = json.loads(OUTPUT.read_text())
        self.assertEqual(data["radii_um"], [150.0, 200.0, 250.0, 300.0])
        self.assertEqual(len(data["patients"]), 8)
        self.assertEqual(len(data["rows"]), 64)
        for result in data["summary_by_radius_um"].values():
            self.assertEqual(result["n_patients"], 8)
            self.assertEqual(len(result["bootstrap_95ci_mean_delta_r2"]), 2)


if __name__ == "__main__":
    unittest.main()
