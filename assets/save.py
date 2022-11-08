'''
Save transformed assembly meshes
'''
import os
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

import numpy as np

from assets.transform import get_transform_matrix


def interpolate_path(path, n_frame=None):
    if n_frame is None: return path
    if len(path) > n_frame:
        interp_path = []
        for i in range(n_frame - 1):
            interp_path.append(path[int(i / n_frame * len(path))])
        interp_path.append(path[-1])
    else:
        interp_path = path
    return interp_path


# def save_path(obj_dir, move_mesh, still_meshes, move_id, still_ids, path, com=np.zeros(3), n_frame=None):
#     '''
#     Save motion of assembly meshes at every time step
#     '''
#     if path is None: return
#     path = interpolate_path(path, n_frame)

#     os.makedirs(obj_dir)
#     for frame, state in enumerate(path):
#         frame_dir = os.path.join(obj_dir, str(frame))
#         os.makedirs(frame_dir)

#         frame_transform = get_transform_matrix(state, com=com)
#         move_mesh_frame = move_mesh.copy()
#         move_mesh_frame.apply_transform(frame_transform)

#         move_obj_path = os.path.join(frame_dir, f'{move_id}.obj')
#         move_mesh_frame.export(move_obj_path, include_color=False, header=None)

#         for still_mesh, still_id in zip(still_meshes, still_ids):
#             still_obj_path = os.path.join(frame_dir, f'{still_id}.obj')
#             still_mesh.export(still_obj_path, include_color=False, header=None)


def save_path(obj_dir, path, com=np.zeros(3), n_frame=None):
    '''
    Save motion of assembly meshes at every time step
    '''
    if path is None: return
    path = interpolate_path(path, n_frame)

    os.makedirs(obj_dir)
    for frame, state in enumerate(path):
        frame_transform = get_transform_matrix(state, com=com)
        np.save(os.path.join(obj_dir, f'{frame}.npy'), frame_transform)


def clear_saved_sdfs(obj_dir):
    for file in os.listdir(obj_dir):
        file_path = os.path.join(obj_dir, file)
        if file_path.endswith('.sdf'):
            os.remove(file_path)
