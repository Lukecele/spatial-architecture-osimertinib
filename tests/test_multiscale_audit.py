import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "MULTISCALE_AUDIT_SUMMARY.json"


class TestMultiscaleAudit(unittest.TestCase):
    def test_summary_preserves_oof_and_corrected_results(self):
        subprocess.run(
            [sys.executable, "src/summarize_multiscale_audit.py"],
            cwd=ROOT,
            check=True,
        )
        data = json.loads(OUTPUT.read_text())
        lung = data["cohorts"]["lung"]
        self.assertEqual(lung["n_scales"], 12)
        self.assertLess(lung["contact_scale_delta_oof_pct"], 0)
        self.assertGreater(lung["best_delta_oof_pct"], 0)
        self.assertTrue(lung["all_in_sample_q_below_0_05"])


if __name__ == "__main__":
    unittest.main()
