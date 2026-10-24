import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_JSON = ROOT / "data" / "public_validation" / "GSE301973_patient_pfs_analysis.json"
OUTPUT_MD = ROOT / "data" / "public_validation" / "GSE301973_PFS_ANALYSIS.md"


class TestGSE301973PatientPFS(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Run or ensure analysis output is fresh
        if not OUTPUT_JSON.exists():
            subprocess.run(
                [sys.executable, "src/analyze_gse301973_patient_pfs.py"],
                cwd=ROOT,
                check=True,
            )

    def test_output_files_exist(self):
        self.assertTrue(OUTPUT_JSON.exists(), f"Missing {OUTPUT_JSON}")
        self.assertTrue(OUTPUT_MD.exists(), f"Missing {OUTPUT_MD}")

    def test_four_patients_mapped_correctly(self):
        data = json.loads(OUTPUT_JSON.read_text())
        self.assertTrue(data.get("mapping_verified"))
        self.assertEqual(len(data["patients"]), 4)

        patient_map = {p["patient_id"]: p for p in data["patients"]}
        expected_ids = {"s1_PFS_69", "s2_PFS_22", "s3_PFS_15", "s4_PFS_12"}
        self.assertEqual(set(patient_map.keys()), expected_ids)

        # Slide A mappings
        self.assertEqual(patient_map["s1_PFS_69"]["slide"], "A")
        self.assertEqual(patient_map["s1_PFS_69"]["sample_position"], "Left")
        self.assertEqual(patient_map["s1_PFS_69"]["PFS_months"], 69.0)
        self.assertEqual(patient_map["s1_PFS_69"]["PFS_group"], "long_PFS")

        self.assertEqual(patient_map["s2_PFS_22"]["slide"], "A")
        self.assertEqual(patient_map["s2_PFS_22"]["sample_position"], "Right")
        self.assertEqual(patient_map["s2_PFS_22"]["PFS_months"], 22.0)
        self.assertEqual(patient_map["s2_PFS_22"]["PFS_group"], "long_PFS")

        # Slide C mappings
        self.assertEqual(patient_map["s3_PFS_15"]["slide"], "C")
        self.assertEqual(patient_map["s3_PFS_15"]["sample_position"], "Left")
        self.assertEqual(patient_map["s3_PFS_15"]["PFS_months"], 15.0)
        self.assertEqual(patient_map["s3_PFS_15"]["PFS_group"], "short_PFS")

        self.assertEqual(patient_map["s4_PFS_12"]["slide"], "C")
        self.assertEqual(patient_map["s4_PFS_12"]["sample_position"], "Right")
        self.assertEqual(patient_map["s4_PFS_12"]["PFS_months"], 12.0)
        self.assertEqual(patient_map["s4_PFS_12"]["PFS_group"], "short_PFS")

    def test_cell_counts_and_confluence(self):
        data = json.loads(OUTPUT_JSON.read_text())
        patients = {p["patient_id"]: p for p in data["patients"]}

        # All patients should have > 10,000 cells
        for pid, p in patients.items():
            self.assertGreater(p["total_cells"], 10000)
            self.assertGreater(p["tumor_cells"], 1000)
            self.assertGreater(p["tumor_fraction"], 0.20)

        # Patient s3_PFS_15 has the highest confluence and homotypic tumor contact
        s3 = patients["s3_PFS_15"]
        self.assertGreater(s3["tumor_fraction"], 0.70)
        self.assertGreater(s3["key_tumor_contacts"]["homotypic_tumor"], 0.80)

        # Long PFS group has lower mean homotypic tumor contact than short PFS
        group_comp = data["group_comparison"]
        mean_homo = group_comp["mean_tumor_homotypic_contact"]
        self.assertLess(mean_homo["long_PFS"], mean_homo["short_PFS"])


if __name__ == "__main__":
    unittest.main()
