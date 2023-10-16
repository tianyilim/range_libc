import numpy as np
import rangelib
import os

'''
This file tests that rangelib can be imported and is able to run.
'''

if __name__ == "__main__":
    FILE_DIR = os.path.dirname(os.path.realpath(__file__))
    TEST_FILE_PATH = os.path.join(
        FILE_DIR, "../../maps/tests/f.png")

    print(TEST_FILE_PATH)

    testOMap = rangelib.OMap(TEST_FILE_PATH)
    testDistanceTransform = rangelib.DistanceTransform(testOMap)
    testRayMarching = rangelib.RayMarching(testDistanceTransform, 50.0)

    query = np.array([[0.1, 0.1, 1.0]], dtype=np.single).transpose()
    print(testRayMarching.batchCalcRange(query))    # should return [50. ]
