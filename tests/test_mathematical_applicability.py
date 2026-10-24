import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from mathematical_applicability_audit import graph_diagnostics


class TestMathematicalApplicability(unittest.TestCase):
    def test_spectral_diagnostics_are_finite(self):
        grid = np.asarray([(x, y) for x in range(5) for y in range(5)], dtype=float)
        signal = grid[:, 0] / 4.0
        result = graph_diagnostics(grid, signal, k=4, n_eigenvalues=8)
        self.assertGreater(result["n_edges"], 0)
        self.assertTrue(np.isfinite(result["truncated_spectral_zeta"]["1"]))
        self.assertTrue(np.isfinite(result["phenotype_dirichlet_energy"]))


if __name__ == "__main__":
    unittest.main()
