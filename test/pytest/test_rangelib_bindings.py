import numpy as np
import rangelib
import os
import unittest

'''
This file tests that rangelib can be imported and is able to run.
'''

FILE_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_FILE_PATH = os.path.join(
    FILE_DIR, "../../maps/tests/f.png")


class TestSum(unittest.TestCase):

    def test_import_rangelib(self):
        """
        Test that rangelib can run
        """

        testOMap = rangelib.OMap(TEST_FILE_PATH)
        testDistanceTransform = rangelib.DistanceTransform(testOMap)
        testRayMarching = rangelib.RayMarching(testDistanceTransform, 50.0)

        query = np.array([[0.1, 0.1, 1.0]], dtype=np.single).transpose()
        res = testRayMarching.batchCalcRange(query)

        self.assertEqual(res.shape[0], 1)
        self.assertAlmostEqual(
            res[0], 50.0
        )


if __name__ == "__main__":
    unittest.main()
