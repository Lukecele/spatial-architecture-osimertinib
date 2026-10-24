import unittest

import numpy as np

from src.rigorous_stats import (
    benjamini_hochberg,
    boolean_mask,
    log_cpm,
    normal_normal_posterior,
    zscore_matrix,
)


class TestRigorousStats(unittest.TestCase):
    def test_boolean_filter_is_explicit(self):
        mask = boolean_mask(
            total_counts=np.array([100, 1000, 2000]),
            detected_features=np.array([300, 100, 500]),
            tissue=np.array([True, True, False]),
        )
        np.testing.assert_array_equal(mask, np.array([False, False, False]))

    def test_matrix_transforms(self):
        matrix = np.array([[10, 20], [10, 20]], dtype=float)
        self.assertEqual(log_cpm(matrix).shape, matrix.shape)
        np.testing.assert_allclose(zscore_matrix(matrix), np.zeros_like(matrix))

    def test_multiple_testing_and_bayes(self):
        adjusted = benjamini_hochberg(np.array([0.01, 0.04, 0.5]))
        self.assertTrue(np.all((adjusted >= 0) & (adjusted <= 1)))
        posterior = normal_normal_posterior(np.array([0.2, 0.3]), 0.1)
        self.assertGreater(posterior["posterior_probability_positive"], 0.5)


if __name__ == "__main__":
    unittest.main()
