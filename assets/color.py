'''
Predefined colors for assembly meshes
'''

import numpy as np


def get_joint_color(index, normalize=True):
    '''
    Get color for 2-part assembly
    '''
    index = int(index)
    assert index in [0, 1]
    colors = np.array([
        [107, 166, 161, 255],
        [209, 184, 148, 255],
    ], dtype=int)
    if normalize: colors = colors.astype(float) / 255.0
    return colors[int(index)]


def get_multi_color(index, normalize=True):
    '''
    Get color for multi-part assembly
    '''
    index = int(index)
    colors = np.array([
        [210, 87, 89, 255],
        [237, 204, 73, 255],
        [60, 167, 221, 255],
        [190, 126, 208, 255],
        [108, 192, 90, 255],
    ], dtype=int)
    if normalize: colors = colors.astype(float) / 255.0
    return colors[index % 5]
