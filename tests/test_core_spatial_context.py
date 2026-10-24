import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "CORE_SPATIAL_CONTEXT_AUDIT.json"


class TestCoreSpatialContext(unittest.TestCase):
    def test_core_audit_has_non_game_models(self):
        subprocess.run(
            [sys.executable, "src/run_core_spatial_context_audit.py"],
            cwd=ROOT,
            check=True,
        )
        data = json.loads(OUTPUT.read_text())
        self.assertEqual(data["radius_um"], 150.0)
        self.assertEqual(set(data["cohorts"]), {"breast", "melanoma", "lung"})
        for result in data["cohorts"].values():
            self.assertIn("focal_composition", result)
            self.assertIn("focal_plus_neighbourhood", result)
            self.assertTrue(-1.0 <= result["focal_plus_neighbourhood"]["spatial_block_oof_r2"] <= 1.0)
            null = result["neighbourhood_permutation_null"]
            self.assertEqual(null["permutations"], 100)
            self.assertGreaterEqual(null["empirical_p_value"], 1 / 101)


if __name__ == "__main__":
    unittest.main()
