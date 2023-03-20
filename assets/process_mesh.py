import os
os.environ['OMP_NUM_THREADS'] = '1'
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

import json
import numpy as np
import trimesh

from assets.subdivide import subdivide_to_size


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


def normalize(meshes, bbox_size=10):

    # compute normalization factors
    vertex_all_stacked = np.vstack([mesh.vertices for mesh in meshes])
    min_box = vertex_all_stacked.min(axis=0)
    max_box = vertex_all_stacked.max(axis=0)
    center = (min_box + max_box) / 2
    scale = max_box - min_box
    scale_factor = bbox_size / np.max(scale)
    scale_transform = norm_to_transform(center, scale_factor)

    # normalization
    for mesh in meshes:
        mesh.apply_transform(scale_transform)
    return meshes


def get_oriented_bounding_box(mesh):
    _, extents = trimesh.bounds.oriented_bounds(mesh)
    return extents


def process_mesh(source_dir, target_dir, subdivide, max_edge=0.5, keep_id=False, verbose=False):
    '''
    1. Check watertight
    2. Scale to unit bounding box
    3. Subdivide (optional)
    4. Translate
    '''
    # order objs
    source_files = []
    for source_file in os.listdir(source_dir):
        if not source_file.endswith('.obj') or source_file == 'assembly.obj': continue
        source_files.append(source_file)
    source_files.sort()

    # load meshes
    meshes = []
    obj_ids = {}
    watertight = True
    for i, source_file in enumerate(source_files):
        if keep_id:
            source_id = int(source_file.replace('.obj', ''))
            obj_ids[source_id] = source_file
        else:
            obj_ids[i] = source_file
        source_path = os.path.join(source_dir, source_file)
        mesh = trimesh.load_mesh(source_path, process=False, maintain_order=True)
        if not mesh.is_watertight:
            print(f'Mesh {source_path} is not watertight')
            watertight = False
        meshes.append(mesh)
    if not watertight:
        return False

    # scale
    meshes = normalize(meshes)

    # subdivide mesh
    if subdivide:
        for i in range(len(meshes)):
            meshes[i] = subdivide_to_size(meshes[i], max_edge=max_edge)

    # compute coms
    coms = {}
    for i, mesh in zip(obj_ids.keys(), meshes):
        com = mesh.center_mass
        coms[i] = com.tolist()
        mesh.apply_transform(com_to_transform(-com))

    # make target directory
    os.makedirs(target_dir, exist_ok=True)

    # write processed objs
    os.makedirs(target_dir, exist_ok=True)
    assert len(meshes) == len(obj_ids)
    for mesh, obj_id in zip(meshes, obj_ids.keys()):
        obj_target_path = os.path.join(target_dir, f'{obj_id}.obj')
        mesh.export(obj_target_path, header=None, include_color=False)
        if verbose:
            print(f'Processed obj written to {obj_target_path}')

    # save translation
    with open(os.path.join(target_dir, 'translation.json'), 'w') as fp:
        json.dump(coms, fp)

    # save obj id map
    with open(os.path.join(target_dir, 'id_map.json'), 'w') as fp:
        json.dump(obj_ids, fp)

    return True


if __name__ == '__main__':
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument('--source-dir', type=str, required=True)
    parser.add_argument('--target-dir', type=str, required=True)
    parser.add_argument('--subdivide', default=False, action='store_true')
    parser.add_argument('--max-edge', type=float, default=0.5)
    parser.add_argument('--keep-id', default=False, action='store_true')
    args = parser.parse_args()

    success = process_mesh(args.source_dir, args.target_dir, args.subdivide, max_edge=args.max_edge, keep_id=args.keep_id, verbose=True)
    print(f'Success: {success}')
