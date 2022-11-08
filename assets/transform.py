'''
Transformation utilities
'''

import numpy as np
from scipy.spatial.transform import Rotation


def get_transform_matrix(state, com=np.zeros(3)):
    '''
    Get transformation matrix of the given state and center of mass
    '''
    if len(state) == 3: # translation only
        transform = np.eye(4)
        transform[:3, 3] = state
        return transform
    elif len(state) == 6: # translation + rotation
        translation, rotation = state[:3], state[3:]
        rotation = Rotation.from_rotvec(rotation).as_matrix()
        trans0_mat = np.eye(4)
        trans0_mat[:3, 3] = -com
        rot_mat = np.eye(4)
        rot_mat[:3, :3] = rotation
        trans1_mat = np.eye(4)
        trans1_mat[:3, 3] = translation + com
        return trans1_mat.dot(rot_mat).dot(trans0_mat)
    else:
        raise NotImplementedError


def get_state_from_matrix(matrix, com=np.zeros(3), full_dof=False):
    '''
    Get state from the given transformation matrix and center of mass
    '''
    translation = matrix[:3, 3]
    if (matrix[:3, :3] == np.eye(3)).all():
        if not full_dof:
            return translation
        else:
            return np.concatenate([translation, np.zeros(3)])
    trans0_mat = np.eye(4)
    trans0_mat[:3, 3] = -com
    rot_mat = np.eye(4)
    rot_mat[:3, :3] = matrix[:3, :3]
    trans1_mat = matrix.dot(np.linalg.inv(rot_mat.dot(trans0_mat)))
    translation = trans1_mat[:3, 3] - com
    rotation = Rotation.from_matrix(rot_mat[:3, :3]).as_rotvec()
    state = np.concatenate([translation, rotation])
    return state


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


def transform_pts_by_state(pts, state, com=np.zeros(3)):
    matrix = get_transform_matrix(state, com)
    return transform_pts_by_matrix(pts, matrix)
