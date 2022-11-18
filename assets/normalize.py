import os
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

import numpy as np
import json
import trimesh
from argparse import ArgumentParser

from assets.load import load_assembly, load_translation


def norm_to_transform(center, scale):
    transform0 = np.eye(4)
    transform0[:3, 3] = -center
    transform1 = np.eye(4)
    transform1[:3, :3] = np.eye(3) * scale
    return transform1.dot(transform0)


def com_to_transform(com):
    transform = np.eye(4)
    transform[:3, 3] = com
    return transform


def normalize(meshes):
    # compute normalization factors
    vertex_all_stacked = np.vstack([mesh.vertices for mesh in meshes])
    min_box = vertex_all_stacked.min(axis=0)
    max_box = vertex_all_stacked.max(axis=0)
    center = (min_box + max_box) / 2
    scale = max_box - min_box
    scale_factor = 10 / np.max(scale)
    scale_transform = norm_to_transform(center, scale_factor)

    # normalization
    for mesh in meshes:
        mesh.apply_transform(scale_transform)
    return meshes


def normalize_multi_assembly(source_dir, target_dir, verbose=False):

    # load meshes
    meshes, names = load_assembly(source_dir, translate=True, return_names=True)
    coms = load_translation(source_dir)

    # normalize
    meshes = normalize(meshes)

    # translation
    coms = {}
    for mesh, name in zip(meshes, names):
        com = mesh.center_mass
        obj_id = int(name.replace('.obj', ''))
        coms[obj_id] = com.tolist()
        mesh.apply_transform(com_to_transform(-com))

    # make target directory
    os.makedirs(target_dir, exist_ok=True)

    # write normalized objs
    os.makedirs(target_dir, exist_ok=True)
    assert len(meshes) == len(names)
    for mesh, name in zip(meshes, names):
        obj_target_path = os.path.join(target_dir, name)
        mesh.export(obj_target_path, header=None, include_color=False)

        if verbose:
            print(f'Normalized obj written to {obj_target_path}')

    # save translation
    with open(os.path.join(target_dir, 'translation.json'), 'w') as fp:
        json.dump(coms, fp)


if __name__ == '__main__':
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument('--source-dir', type=str, required=True)
    parser.add_argument('--target-dir', type=str, required=True)
    args = parser.parse_args()

    normalize_multi_assembly(args.source_dir, args.target_dir, verbose=True)
