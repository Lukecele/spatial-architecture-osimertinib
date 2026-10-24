import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "public_validation" / "GSE176078_validation.json"


class TestGSE176078Validation(unittest.TestCase):
    def test_validation_output_is_patient_level(self):
        data = json.loads(OUTPUT.read_text())
        self.assertEqual(data["accession"], "GSE176078")
        self.assertEqual(len(data["patients"]), 6)
        self.assertEqual(data["radius_um"], 150.0)
        self.assertEqual(len(data["bootstrap_95ci_mean_delta_r2"]), 2)
        self.assertTrue(all("context_minus_focal_r2" in row for row in data["patients"]))


if __name__ == "__main__":
    unittest.main()
