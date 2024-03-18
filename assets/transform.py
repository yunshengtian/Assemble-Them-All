'''
Transformation utilities
'''

import numpy as np
from scipy.spatial.transform import Rotation


def get_transform_matrix(state):
    '''
    Get transformation matrix of the given state
    '''
    if len(state) == 3: # translation only
        transform = np.eye(4)
        transform[:3, 3] = state
        return transform
    elif len(state) == 6: # translation + rotation
        translation, rotation = state[:3], state[3:]
        rotation = Rotation.from_rotvec(rotation).as_matrix()
        transform = np.eye(4)
        transform[:3, :3] = rotation
        transform[:3, 3] = translation
        return transform
    else:
        raise NotImplementedError


def transform_pts_by_matrix(pts, matrix):
    '''
    Transform an array of xyz pts (n, 3) by a 4x4 matrix
    '''
    pts = np.array(pts)
    if len(pts.shape) == 1:
        if len(pts) == 3:
            v = np.append(pts, 1.0)
        elif len(pts) == 4:
            v = pts
        else:
            raise NotImplementedError
        v = matrix @ v
        return v[0:3]
    elif len(pts.shape) == 2:
        # transpose first
        if pts.shape[1] == 3:
            # pad the points with ones to be (n, 4)
            v = np.hstack([pts, np.ones((len(pts), 1))]).T
        elif pts.shape[1] == 4:
            v = pts.T
        else:
            raise NotImplementedError
        v = matrix @ v
        # transpose and crop back to (n, 3)
        return v.T[:, 0:3]
    else:
        raise NotImplementedError


def transform_pts_by_state(pts, state):
    matrix = get_transform_matrix(state)
    return transform_pts_by_matrix(pts, matrix)
