import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "public_validation" / "GSE189487_validation.json"


class TestGSE189487Validation(unittest.TestCase):
    def test_lung_validation_output(self):
        data = json.loads(OUTPUT.read_text())
        self.assertEqual(data["accession"], "GSE189487")
        self.assertEqual(len(data["patients"]), 6)
        self.assertEqual(data["radius_um"], 150.0)
        self.assertEqual(len(data["bootstrap_95ci_mean_delta_r2"]), 2)
        self.assertEqual(len(data["marker_spec"]["immune"]), 5)


if __name__ == "__main__":
    unittest.main()
