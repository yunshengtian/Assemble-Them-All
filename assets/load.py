'''
Load assembly meshes and translations
'''
import os
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

import numpy as np
import json
import trimesh

from assets.color import get_joint_color, get_multi_color
from assets.transform import get_transform_matrix, transform_pts_by_state


def load_translation(obj_dir, rotvec=None):
    '''
    Load translation from dir
    '''
    coms = None
    translation_path = os.path.join(obj_dir, 'translation.json')
    if rotvec is not None:
        assert len(rotvec) == 3
    if os.path.exists(translation_path):
        with open(translation_path, 'r') as fp:
            coms = json.load(fp)
        new_coms = {}
        for key, val in coms.items():
            new_coms[int(key)] = np.array(val)
            if rotvec is not None:
                new_coms[int(key)] = transform_pts_by_state(new_coms[int(key)], np.concatenate([np.zeros(3), rotvec]))
        coms = new_coms
    return coms


def com_to_transform(com):
    '''
    COM to transformation matrix
    '''
    transform = np.eye(4)
    transform[:3, 3] = com
    return transform


def load_assembly(obj_dir, translate=True, rotvec=None, return_names=False):
    '''
    Load the entire assembly from dir
    '''
    meshes = []
    names = []
    if translate:
        coms = load_translation(obj_dir)
    else:
        coms = None

    obj_ids = []
    for file_name in os.listdir(obj_dir):
        if file_name.endswith('.obj'):
            try:
                obj_id = int(file_name.replace('.obj', ''))
            except:
                continue
            obj_ids.append(obj_id)
    obj_ids = sorted(obj_ids)
    get_color = get_joint_color if len(obj_ids) <= 2 else get_multi_color

    for seq, obj_id in enumerate(obj_ids):
        obj_name = f'{obj_id}.obj'
        obj_path = os.path.join(obj_dir, obj_name)
        mesh = trimesh.load_mesh(obj_path, process=False, maintain_order=True)
        if rotvec is not None:
            assert len(rotvec) == 3
            rot_transform = get_transform_matrix(np.concatenate([np.zeros(3), rotvec]))
            mesh.apply_transform(rot_transform)
        if coms is not None:
            if rotvec is not None:
                coms[obj_id] = transform_pts_by_state(coms[obj_id], np.concatenate([np.zeros(3), rotvec]))
            com_transform = com_to_transform(coms[obj_id])
            mesh.apply_transform(com_transform)
        mesh.visual.face_colors = get_color(seq, normalize=False)
        meshes.append(mesh)
        names.append(obj_name)
    
    if return_names:
        return meshes, names
    else:
        return meshes


def load_paths(path_dir):
    '''
    Load motion of assembly meshes at every time step
    '''
    paths = {}
    for step in os.listdir(path_dir):
        obj_id = int(step.split('_')[1])
        step_dir = os.path.join(path_dir, step)
        if os.path.isdir(step_dir):
            path = []
            frame_files = []
            for frame_file in os.listdir(step_dir):
                if frame_file.endswith('.npy'):
                    frame_files.append(frame_file)
            frame_files.sort(key=lambda x: int(x.replace('.npy', '')))
            for frame_file in frame_files:
                frame_path = os.path.join(step_dir, frame_file)
                frame_transform = np.load(frame_path)
                path.append(frame_transform)
            paths[obj_id] = path
    return paths


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--dir', type=str, required=True)
    args = parser.parse_args()

    meshes = load_assembly(args.dir)
    trimesh.Scene(meshes).show()
