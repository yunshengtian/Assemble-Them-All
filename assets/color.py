'''
Predefined colors for assembly meshes
'''

import numpy as np


def get_color(part_ids, normalize=True):
    color_map = {}
    if len(part_ids) <= 2:
        colors = np.array([
            [107, 166, 161, 255],
            [209, 184, 148, 255],
        ], dtype=int)
    else:
        colors = np.array([
            [210, 87, 89, 255],
            [237, 204, 73, 255],
            [60, 167, 221, 255],
            [190, 126, 208, 255],
            [108, 192, 90, 255],
        ], dtype=int)
    if normalize: colors = colors.astype(float) / 255.0
    for i, part_id in enumerate(part_ids):
        color_map[part_id] = colors[i % len(colors)]
    return color_map
