'''
Subdivide a mesh until every edge is shorter than a specified length
'''

import os
import shutil
import numpy as np
import trimesh
from sortedcontainers import SortedKeyList

from load import load_assembly


def subdivide_to_size(mesh, max_edge):
    vertices = mesh.vertices
    faces = mesh.faces
    edges = mesh.edges

    # check if reversed edges exist
    edges_rev = edges[:, ::-1]
    edges_unique = np.unique(np.vstack([edges, edges_rev]), axis=0)
    assert len(edges_unique) == len(edges)

    # quick reference from edge to face and length
    edges_face = mesh.edges_face
    edges_length = np.linalg.norm(vertices[edges][:, 0] - vertices[edges][:, 1], axis=1)
    edge_face_dict = dict(zip(map(tuple, edges), edges_face))
    edge_length_dict = dict(zip(map(tuple, edges), edges_length))
    edge_sorted = SortedKeyList(map(tuple, edges), key=lambda e: edge_length_dict[e])

    while True:
        # check longest edge
        longest_edge = edge_sorted[-1]
        longest_edge_length = edge_length_dict[longest_edge]
        if longest_edge_length < max_edge:
            break

        # build 1 new vertex
        new_vertex = vertices[list(longest_edge)].mean(axis=0)

        # add 1 new vertex
        new_vertex_idx = len(vertices)
        vertices = np.vstack([vertices, new_vertex])

        # build 4 new faces
        vertex0_idx, vertex1_idx = longest_edge
        face0_idx, face1_idx = edge_face_dict[longest_edge], edge_face_dict[longest_edge[::-1]]
        face0, face1 = faces[face0_idx], faces[face1_idx]
        new_face00, new_face01 = face0.copy(), face0.copy()
        new_face10, new_face11 = face1.copy(), face1.copy()
        new_face00[face0 == vertex0_idx] = new_vertex_idx
        new_face01[face0 == vertex1_idx] = new_vertex_idx
        new_face10[face1 == vertex0_idx] = new_vertex_idx
        new_face11[face1 == vertex1_idx] = new_vertex_idx

        # remove 2 old faces and add 4 new faces
        new_face00_idx, new_face01_idx = face0_idx, face1_idx
        new_face10_idx, new_face11_idx = len(faces), len(faces) + 1
        faces[face0_idx] = new_face00
        faces[face1_idx] = new_face01
        faces = np.vstack([faces, [new_face10, new_face11]])

        # remove 2 old edges and add 8 new edges
        edge_sorted.remove(longest_edge)
        edge_sorted.remove(longest_edge[::-1])
        edge_face_dict.pop(longest_edge)
        edge_face_dict.pop(longest_edge[::-1])
        edge_length_dict.pop(longest_edge)
        edge_length_dict.pop(longest_edge[::-1])

        for new_face, new_face_idx in zip([new_face00, new_face01, new_face10, new_face11], [new_face00_idx, new_face01_idx, new_face10_idx, new_face11_idx]):
            new_edges = [tuple(new_face[[0, 1]]), tuple(new_face[[1, 2]]), tuple(new_face[[2, 0]])]
            for new_edge in new_edges:
                edge_face_dict[new_edge] = new_face_idx
                edge_length_dict[new_edge] = np.linalg.norm(vertices[new_edge[0]] - vertices[new_edge[1]])
                edge_sorted.discard(new_edge)
                edge_sorted.add(new_edge)

    result = trimesh.Trimesh(
            vertices=vertices,
            faces=faces,
            process=False,
            maintain_order=True)
    assert result.is_watertight
    return result


def subdivide_assembly(source_dir, target_dir, max_edge, render=False):
    # load objs
    meshes, names = load_assembly(source_dir, translate=False, return_names=True)

    if render:
        trimesh.Scene(meshes).show()

    # subdivide
    new_meshes = []
    for mesh in meshes:
        new_mesh = subdivide_to_size(mesh, max_edge)
        new_meshes.append(new_mesh)

    if render:
        trimesh.Scene(new_meshes).show()

    # export objs
    os.makedirs(target_dir, exist_ok=True)
    for new_mesh, name in zip(new_meshes, names):
        target_obj_path = os.path.join(target_dir, name)
        new_mesh.export(target_obj_path, header=None)

    # copy translation
    source_trans_path = os.path.join(source_dir, 'translation.json')
    target_trans_path = os.path.join(target_dir, 'translation.json')
    if os.path.exists(source_trans_path):
        shutil.copyfile(source_trans_path, target_trans_path)


if __name__ == '__main__':
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument('--source-dir', type=str, required=True)
    parser.add_argument('--target-dir', type=str, required=True)
    parser.add_argument('--max-edge', type=float, default=0.5)
    parser.add_argument('--render', default=False, action='store_true')
    args = parser.parse_args()

    subdivide_assembly(args.source_dir, args.target_dir, args.max_edge, render=args.render)
