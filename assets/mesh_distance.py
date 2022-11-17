'''
Compute the minimum distance between meshes at certain states
'''

import os
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

import numpy as np
from assets.transform import get_transform_matrix, transform_pts_by_matrix


def compute_all_mesh_distance(meshes, states, coms=None):
    '''
    Compute the minimum distance between meshes at certain states
    '''
    assert len(meshes) == len(states)
    if coms is None:
        coms = [np.zeros(3) for _ in meshes]
    else:
        assert len(coms) == len(meshes)

    mats, inv_mats = [], []
    for i in range(len(meshes)):
        mat = get_transform_matrix(states[i], coms[i])
        mats.append(mat)
        inv_mats.append(np.linalg.inv(mat))

    d = np.inf
    for i in range(len(meshes)):
        for j in range(i + 1, len(meshes)):
            v_i_trans = transform_pts_by_matrix(meshes[i].vertices.T, inv_mats[j].dot(mats[i]))
            v_j_trans = transform_pts_by_matrix(meshes[j].vertices.T, inv_mats[i].dot(mats[j]))
            d_ij = meshes[i].min_distance(v_j_trans.T)
            d_ji = meshes[j].min_distance(v_i_trans.T)
            d = np.min([d, d_ij, d_ji])
    
    return d


def compute_move_mesh_distance(move_mesh, still_meshes, state, com=None):
    '''
    Compute the minimum distance between meshes at certain states
    '''
    if com is None: com = np.zeros(3)

    move_mat = get_transform_matrix(state, com)
    move_inv_mat = np.linalg.inv(move_mat)

    v_m_trans = transform_pts_by_matrix(move_mesh.vertices.T, move_mat).T

    d = np.inf
    for i in range(len(still_meshes)):
        v_s_trans = transform_pts_by_matrix(still_meshes[i].vertices.T, move_inv_mat).T
        d_ms = move_mesh.min_distance(v_s_trans)
        d_sm = still_meshes[i].min_distance(v_m_trans)
        d = np.min([d, d_ms, d_sm])

    return d
