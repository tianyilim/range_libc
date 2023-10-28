'''
This file generates some of the maps used later for benchmarking.
'''

import os
import numpy as np
from PIL import Image

FILE_DIR = os.path.dirname(os.path.realpath(__file__))
MAPS_FILE_PATH = os.path.join(
    FILE_DIR, "../maps/tests")

MAP_SIZE = 600  # we want 600x600 square maps

corridor_map = np.zeros((MAP_SIZE, MAP_SIZE), dtype=np.uint8)
corridor_map[150:450, :] = 255
Image.fromarray(corridor_map).save(
    os.path.join(MAPS_FILE_PATH, "corridor.png"))

corner_map = np.zeros((MAP_SIZE, MAP_SIZE), dtype=np.uint8)
corner_map[150:450, 150:450] = 255
corner_map[150:450, 0:150] = 255
corner_map[450:, 150:450] = 255
Image.fromarray(corner_map).save(
    os.path.join(MAPS_FILE_PATH, "corner.png"))

# curve_map = np.zeros((MAP_SIZE, MAP_SIZE), dtype=np.uint8)
