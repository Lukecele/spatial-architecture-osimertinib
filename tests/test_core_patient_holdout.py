import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "extension_pilot" / "core_patient_holdout.json"


class TestCorePatientHoldout(unittest.TestCase):
    def test_core_patient_holdout_has_eight_patients(self):
        subprocess.run(
            [sys.executable, "src/core_patient_holdout.py"],
            cwd=ROOT,
            check=True,
        )
        data = json.loads(OUTPUT.read_text())
        self.assertEqual(data["radius_um"], 150.0)
        self.assertEqual(len(data["patients"]), 8)
        self.assertEqual(len(data["patient_summary"]), 8)
        self.assertEqual(len(data["bootstrap_95ci_mean_delta_r2"]), 2)


if __name__ == "__main__":
    unittest.main()
