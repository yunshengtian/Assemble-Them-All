import os
os.environ['OMP_NUM_THREADS'] = '1'
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

import numpy as np
import redmax_py as redmax
from argparse import ArgumentParser
import trimesh
from time import time
import random
from scipy.spatial.transform import Rotation

from pyplanners.rrt_connect import rrt_connect, birrt
from pyplanners.rrt import rrt
from pyplanners.targetless_rrt import targetless_rrt
from pyplanners.smoothing import smooth_path

from assets.load import load_assembly, load_part_ids
from assets.save import interpolate_path, save_path
from assets.color import get_color
from assets.transform import transform_pts_by_state, get_transform_matrix
from assets.mesh_distance import compute_move_mesh_distance


class PyPlanner:

    def __init__(self, assembly_dir, move_id, still_ids, rotation=False, max_collision=1e-2, adaptive_collision=False, body_type='bvh', sdf_dx=0.05, save_sdf=False):
        '''
        Initialize meshes and their geometric properties
        '''
        # params
        self.rotation = rotation
        self.max_collision = max_collision

        # load data
        meshes, names = load_assembly(assembly_dir, return_names=True)
        self.move_id, self.still_ids = move_id, still_ids
        part_ids = load_part_ids(assembly_dir)
        color_map = get_color(part_ids, normalize=False)

        # build meshes
        self.viz_mesh_move, self.viz_meshes_still = None, []
        self.mesh_move, self.meshes_still = None, []

        for mesh, name in zip(meshes, names):
            mesh_id = name.replace('.obj', '')
            mesh.visual.face_colors = color_map[mesh_id]

            sdf_load_path = os.path.join(assembly_dir, name.replace('.obj', '.sdf'))
            sdf_save_path = sdf_load_path if save_sdf else ''

            if mesh_id == move_id:
                self.viz_mesh_move = mesh
                if body_type == 'bvh':
                    self.mesh_move = redmax.BVHMesh(mesh.vertices.T, mesh.faces.T)
                elif body_type == 'sdf':
                    self.mesh_move = redmax.SDFMesh(mesh.vertices.T, mesh.faces.T, sdf_dx, sdf_load_path, sdf_save_path)
                else:
                    raise NotImplementedError

            elif mesh_id in still_ids:
                self.viz_meshes_still.append(mesh)
                if body_type == 'bvh':
                    self.meshes_still.append(redmax.BVHMesh(mesh.vertices.T, mesh.faces.T))
                elif body_type == 'sdf':
                    self.meshes_still.append(redmax.SDFMesh(mesh.vertices.T, mesh.faces.T, sdf_dx, sdf_load_path, sdf_save_path))
                else:
                    raise NotImplementedError

        if adaptive_collision:
            min_d = compute_move_mesh_distance(self.mesh_move, self.meshes_still, state=np.zeros(3))
            self.max_collision = max(-min_d, 0) + self.max_collision

        # compute bounding boxes
        self.vertices_move = self.viz_mesh_move.vertices
        self.min_box_move = self.vertices_move.min(axis=0)
        self.max_box_move = self.vertices_move.max(axis=0)
        self.size_box_move = self.max_box_move - self.min_box_move

        self.vertices_still = np.vstack([mesh.vertices for mesh in self.viz_meshes_still])
        self.min_box_still = self.vertices_still.min(axis=0)
        self.max_box_still = self.vertices_still.max(axis=0)
        self.size_box_still = self.max_box_still - self.min_box_still

        # compute bounds
        self.state_lower_bound = (self.min_box_still - self.max_box_move) - 0.5 * self.size_box_move
        self.state_upper_bound = (self.max_box_still - self.min_box_move) + 0.5 * self.size_box_move

        # compute convex hulls
        self.hull_move = trimesh.convex.convex_hull(self.vertices_move)
        self.hull_still = trimesh.convex.convex_hull(self.vertices_still)
        self.collision_manager = trimesh.collision.CollisionManager()
        self.collision_manager.add_object('hull_still', self.hull_still)

    def get_distance_fn(self):
        if self.rotation:
            def distance_fn(q1, q2):
                boxes1 = transform_pts_by_state(np.vstack([self.min_box_move, self.max_box_move]), q1)
                boxes2 = transform_pts_by_state(np.vstack([self.min_box_move, self.max_box_move]), q2)
                return np.linalg.norm(boxes1 - boxes2, axis=1).sum()
        else:
            def distance_fn(q1, q2):
                return np.linalg.norm(q1 - q2)
        return distance_fn

    def get_collision_fn(self):
        def collision_fn(q):
            d = compute_move_mesh_distance(self.mesh_move, self.meshes_still, q)
            return d < -self.max_collision
        return collision_fn

    def get_sample_fn(self, collision_fn):
        if self.rotation:
            def sample_fn():
                while True:
                    ratio = np.random.random(3)
                    translation = ratio * self.state_lower_bound + (1.0 - ratio) * self.state_upper_bound
                    rotation = Rotation.random().as_rotvec()
                    state = np.concatenate([translation, rotation])
                    if not collision_fn(state):
                        return state
        else:
            def sample_fn():
                while True:
                    ratio = np.random.random(3)
                    state = ratio * self.state_lower_bound + (1.0 - ratio) * self.state_upper_bound
                    if not collision_fn(state):
                        return state
        return sample_fn

    def get_extend_fn(self, distance_fn, step_size):
        def extend_fn(q1, q2):
            q_dist = distance_fn(q1, q2)
            num_steps = int(np.ceil(q_dist / step_size))
            for i in range(num_steps):
                yield q1 + (q2 - q1) / q_dist * (i + 1) * step_size
        return extend_fn

    def get_goal_test(self):
        def goal_test(q):
            transform = get_transform_matrix(q)
            hull_move = self.hull_move.copy()
            hull_move.apply_transform(transform)
            has_collision = self.collision_manager.in_collision_single(hull_move)
            if not has_collision:
                min_box_move, max_box_move = hull_move.vertices.min(axis=0), hull_move.vertices.max(axis=0)
                move_contain_still = (min_box_move <= self.min_box_still).all() and (max_box_move >= self.max_box_still).all()
                still_contain_move = (self.min_box_still <= min_box_move).all() and (self.max_box_still >= max_box_move).all()
                return not (move_contain_still or still_contain_move) # check if one hull fully contains another
            else:
                return False
        return goal_test

    def visualize_state(self, state, rotvec=[0, 0, 0]):
        viz_mesh_move = self.viz_mesh_move.copy()
        viz_meshes_still = [viz_mesh_still.copy() for viz_mesh_still in self.viz_meshes_still]
        viz_mesh_move.vertices = transform_pts_by_state(viz_mesh_move.vertices, state)

        assert len(rotvec) == 3
        transform = np.eye(4)
        transform[:3, :3] = Rotation.from_rotvec(rotvec).as_matrix()

        viz_mesh_move.apply_transform(transform)
        for viz_mesh_still in viz_meshes_still:
            viz_mesh_still.apply_transform(transform)
        trimesh.Scene([viz_mesh_move, *viz_meshes_still]).show()

    def visualize_path(self, path, n_viz_state=4, rotvec=[0, 0, 0]):
        viz_path = interpolate_path(path, n_viz_state)
        for state in viz_path:
            self.visualize_state(state, rotvec)

    def save_path(self, path, save_dir, n_save_state):
        save_path(save_dir, self.viz_mesh_move, self.viz_meshes_still, self.move_id, self.still_ids, path, n_frame=n_save_state)

    def seed(self, seed):
        random.seed(seed)
        np.random.seed(seed)

    def plan(self, planner_name, step_size, max_time, seed=1, return_path=False, simplify=False, render=False):

        self.seed(seed)

        distance_fn = self.get_distance_fn()
        collision_fn = self.get_collision_fn()
        sample_fn = self.get_sample_fn(collision_fn)
        extend_fn = self.get_extend_fn(distance_fn, step_size)
        goal_test = self.get_goal_test()

        if self.rotation:
            start = np.zeros(6, dtype=float)
        else:
            start = np.zeros(3, dtype=float)

        if collision_fn(start):
            status = 'Invalid start state'
            return (status, 0., None) if return_path else (status, 0.)

        if goal_test(start):
            status = 'Start with goal'
            return (status, 0., None) if return_path else (status, 0.)

        while True:
            goal = sample_fn()
            if goal_test(goal):
                break

        if render:
            self.visualize_state(start)
            self.visualize_state(goal)

        if planner_name == 'matevec-trrt':

            # compute mating vectors
            face_normals = np.vstack([mesh.face_normals for mesh in [self.viz_mesh_move, *self.viz_meshes_still]])
            face_normals = np.unique(face_normals, axis=0)

            mating_vectors = np.empty((0, 3), dtype=float)
            for i in range(len(face_normals)):
                if (mating_vectors @ face_normals[i] > 0.95).any():
                    pass
                else:
                    mating_vectors = np.vstack([mating_vectors, face_normals[i]])

        status = 'Failure'
        t_start = time()

        path = None
        max_iterations = sys.maxsize
        if planner_name == 'rrt':
            path = rrt(start, goal, distance_fn, sample_fn, extend_fn, collision_fn, goal_test=goal_test, max_iterations=max_iterations, max_time=max_time)
        elif planner_name == 'rrt-connect':
            path = rrt_connect(start, goal, distance_fn, sample_fn, extend_fn, collision_fn, max_iterations=max_iterations, max_time=max_time)
        elif planner_name == 'birrt':
            path = birrt(start, goal, distance_fn, sample_fn, extend_fn, collision_fn, max_iterations=max_iterations, max_time=max_time)
        elif planner_name == 'trrt':
            path = targetless_rrt(start, self.vertices_move, self.vertices_still,
                distance_fn, sample_fn, extend_fn, collision_fn, goal_test=goal_test,
                goal_probability=0.2, max_iterations=max_iterations, max_time=max_time)
        elif planner_name == 'matevec-trrt':
            path = None
            if not self.rotation:
                for mating_vector in mating_vectors:
                    for q in extend_fn(start, start + mating_vector * 20):
                        if time() - t_start > max_time:
                            status = 'Timeout'
                            break
                        if collision_fn(q):
                            break
                        if goal_test(q):
                            path = [start, q]
                            status = 'Success'
                            break
                    if status != 'Failure':
                        break
            if status == 'Failure':
                path = targetless_rrt(start, self.vertices_move, self.vertices_still,
                    distance_fn, sample_fn, extend_fn, collision_fn, goal_test=goal_test,
                    goal_probability=0.2, max_iterations=max_iterations, max_time=max_time - (time() - t_start))
        else:
            raise NotImplementedError

        t_plan = time() - t_start

        status = 'Success' if path is not None else 'Timeout'

        if status == 'Success':

            if simplify:
                path = smooth_path(path, extend_fn, collision_fn, distance_fn, None, sample_fn, max_time=max_time)

            if render:
                self.visualize_path(path)

        return (status, t_plan, path) if return_path else (status, t_plan)
    

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--dir', type=str, default='joint_assembly', help='directory storing all assemblies')
    parser.add_argument('--id', type=str, required=True, help='assembly id (e.g. 00000)')
    parser.add_argument('--move-id', type=str, default='0')
    parser.add_argument('--still-ids', type=str, nargs='+', default=['1'])
    parser.add_argument('--rotation', default=False, action='store_true')
    parser.add_argument('--planner', type=str, required=True, choices=['rrt', 'rrt-connect', 'birrt', 'trrt', 'matevec-trrt'])
    parser.add_argument('--max-collision', type=float, default=1e-2)
    parser.add_argument('--adaptive-collision', default=False, action='store_true')
    parser.add_argument('--step-size', type=float, default=1e-2)
    parser.add_argument('--max-time', type=float, default=120)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--body-type', type=str, default='sdf', choices=['bvh', 'sdf'], help='simulation type of body')
    parser.add_argument('--sdf-dx', type=float, default=0.05, help='grid resolution of SDF')
    parser.add_argument('--simplify', default=False, action='store_true')
    parser.add_argument('--render', default=False, action='store_true')
    parser.add_argument('--save-dir', type=str, default=None)
    parser.add_argument('--n-save-state', type=int, default=100)
    args = parser.parse_args()

    asset_folder = os.path.join(project_base_dir, './assets')
    assembly_dir = os.path.join(asset_folder, args.dir, args.id)

    planner = PyPlanner(assembly_dir, args.move_id, args.still_ids, rotation=args.rotation, max_collision=args.max_collision, adaptive_collision=args.adaptive_collision,
        body_type=args.body_type, sdf_dx=args.sdf_dx)
    status, t_plan, path = planner.plan(
        planner_name=args.planner, 
        step_size=args.step_size, 
        max_time=args.max_time, 
        seed=args.seed, 
        return_path=True, 
        simplify=args.simplify, 
        render=args.render
    )

    print(f'Status: {status}, planning time: {t_plan}')

    if args.save_dir is not None:
        planner.save_path(path, args.save_dir, args.n_save_state)