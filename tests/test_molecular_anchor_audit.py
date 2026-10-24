import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "MOLECULAR_ANCHOR_AUDIT.json"


class TestMolecularAnchorAudit(unittest.TestCase):
    def test_missing_mirna_and_methylome_are_explicit(self):
        subprocess.run(
            [sys.executable, "src/molecular_anchor_audit.py"],
            cwd=ROOT,
            check=True,
        )
        data = json.loads(OUTPUT.read_text())
        melanoma = data["anchors"]["melanoma_MITF_APAF1_resistance"]
        self.assertEqual(melanoma["present_features"], ["MITF"])
        self.assertEqual(melanoma["missing_features"], ["APAF1"])
        self.assertIn("not identifiable", data["anchors"]["PDAC_miR1307_CLIC5_FOLFIRINOX"]["status"])


if __name__ == "__main__":
    unittest.main()
