import os
os.environ['OMP_NUM_THREADS'] = '1'
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

import numpy as np
import redmax_py as redmax
import json
import trimesh
from time import time
import networkx as nx
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import pdist, squareform
from pyquaternion import Quaternion

from assets.load import load_translation, load_assembly
from assets.save import save_path, clear_saved_sdfs
from assets.color import get_joint_color, get_multi_color
from assets.transform import transform_pts_by_state
from assets.mesh_distance import compute_move_mesh_distance
from utils.renderer import SimRenderer


class State:

    def __init__(self, q, qdot):
        self.q = q
        self.qdot = qdot

    def __repr__(self):
        return f'[State object at {hex(id(self))}]'


class Tree:

    def __init__(self):
        self.tree = nx.DiGraph()
        self.node_idx = 0
        self.root_node = None
        self.last_node = None

    def add_node(self, state):
        self.tree.add_node(state, idx=self.node_idx)
        self.node_idx += 1
        if self.node_idx == 1:
            self.root_node = state
        self.last_node = state

    def add_edge(self, state_src, state_tar, action, states_between):
        assert self.tree.has_node(state_tar)
        assert self.tree.has_node(state_tar)
        self.tree.add_edge(state_src, state_tar, action=action, states_between=states_between)

    def get_in_edge(self, state):
        assert self.tree.has_node(state)
        in_edges = list(self.tree.in_edges(state, data=True))
        assert len(in_edges) <= 1
        if len(in_edges) == 0:
            return None
        else:
            return in_edges[0]

    def get_out_edges(self, state):
        assert self.tree.has_node(state)
        return list(self.tree.out_edges(state, data=True))

    def has_path(self, start_state, end_state):
        state = end_state
        while True:
            if state == start_state:
                return True
            in_edge = self.get_in_edge(state)
            if in_edge is None:
                return False
            else:
                state = in_edge[0]

    def get_path(self, start_state, end_state):
        path = []
        state = end_state
        while True:
            path.append(state)
            if state == start_state:
                break
            in_edge = self.get_in_edge(state)
            if in_edge is None:
                return None # failed
            else:
                state_pred, states_betweeen = in_edge[0], in_edge[2]['states_between']
                path.extend(states_betweeen[::-1])
                state = state_pred
        return path[::-1]

    def get_root_path(self, state):
        return self.get_path(self.root_node, state)

    def get_nodes(self):
        return list(self.tree.nodes)

    def get_edges(self):
        return list(self.tree.edges(data=True))

    def get_predecessor(self, state):
        preds = list(self.tree._pred[state].keys())
        assert len(preds) <= 1
        if len(preds) == 0:
            return None
        else:
            return preds[0]

    def get_successors(self, state):
        return list(self.tree._succ[state].keys())

    def get_in_degree(self, state):
        return self.tree.in_degree(state)

    def get_out_degree(self, state):
        return self.tree.out_degree(state)

    def get_node_attr(self, node, attr_name):
        return self.tree.nodes[node][attr_name]

    def set_node_attr(self, node, attr_name, attr_val):
        self.tree.nodes[node][attr_name] = attr_val

    def draw(self):
        from networkx.drawing.nx_pydot import graphviz_layout
        import matplotlib.pyplot as plt
        pos = graphviz_layout(self.tree, prog='dot')
        nx.draw_networkx(self.tree, pos, arrows=True, with_labels=False, node_size=10)
        plt.show()


def arr_to_str(arr):
    return ' '.join([str(x) for x in arr])


def get_xml_string(assembly_dir, move_id, still_ids, move_joint_type, body_type, sdf_dx, col_th, save_sdf):
    with open(os.path.join(assembly_dir, 'translation.json'), 'r') as fp:
        translation = json.load(fp)
    body_type = body_type.upper()
    get_color = get_joint_color if len(translation.keys()) <= 2 else get_multi_color
    sdf_args = 'load_sdf="true" save_sdf="true"' if save_sdf else ''
    string = f'''
<redmax model="assemble">
<option integrator="BDF1" timestep="1e-3" gravity="0. 0. 1e-12"/>
<render bg_color="255 255 255 255"/>

<default>
    <general_{body_type}_contact kn="1e6" kt="1e3" mu="0" damping="0"/>
</default>
'''
    for part_id in [move_id, *still_ids]:
        joint_type = move_joint_type if part_id == move_id else 'fixed'
        string += f'''
<robot>
    <link name="part{part_id}">
        <joint name="part{part_id}" type="{joint_type}" axis="0. 0. 0." pos="{arr_to_str(translation[str(part_id)])}" quat="1 0 0 0" frame="WORLD" damping="0"/>
        <body name="part{part_id}" type="{body_type}" filename="{assembly_dir}/{part_id}.obj" {sdf_args} pos="0 0 0" quat="1 0 0 0" scale="1 1 1" transform_type="OBJ_TO_JOINT" density="1" dx="{sdf_dx}" col_th="{col_th}" mu="0" rgba="{arr_to_str(get_color(part_id))}"/>
    </link>
</robot>
'''
    string += f'''
<contact>
'''
    for part_id in still_ids:
        string += f'''
    <general_{body_type}_contact general_body="part{move_id}" {body_type}_body="part{part_id}"/>
    <general_{body_type}_contact general_body="part{part_id}" {body_type}_body="part{move_id}"/>
'''
    string += f'''
</contact>
</redmax>
    '''
    return string


def unit_vector(vector):
    norm = np.linalg.norm(vector)
    if np.isclose(norm, 0):
        random_vector = np.random.random(len(vector))
        return random_vector / np.linalg.norm(random_vector)
    else:
        return vector / norm


class PhysicsPlanner:

    def __init__(self, asset_folder, assembly_dir, move_id, still_ids, 
        rotation=False, body_type='bvh', sdf_dx=0.05, collision_th=0.01, force_mag=1e3, frame_skip=100, save_sdf=False):

        # calculate collision threshold
        meshes, names = load_assembly(assembly_dir, return_names=True)
        move_mesh = None
        still_meshes = []
        for mesh, name in zip(meshes, names):
            if body_type == 'bvh':
                phys_mesh = redmax.BVHMesh(mesh.vertices.T, mesh.faces.T)
            elif body_type == 'sdf':
                phys_mesh = redmax.SDFMesh(mesh.vertices.T, mesh.faces.T, sdf_dx, load_path=os.path.join(assembly_dir, name.replace('.obj', '.sdf')))
            else:
                raise NotImplementedError
            if name == f'{move_id}.obj':
                move_mesh = phys_mesh
            else:
                still_meshes.append(phys_mesh)
        min_d = compute_move_mesh_distance(move_mesh, still_meshes, state=np.zeros(3))
        collision_th = max(-min_d, 0) + collision_th
        
        # build simulation
        move_joint_type = 'free3d-exp' if rotation else 'translational'
        model_string = get_xml_string(assembly_dir, move_id, still_ids, move_joint_type, body_type, sdf_dx, collision_th, save_sdf)
        self.sim = redmax.Simulation(model_string, asset_folder)
        self.rotation = rotation
        self.ndof = 6 if rotation else 3
        self.force_mag = force_mag
        self.frame_skip = frame_skip

        # names
        self.move_id, self.still_ids = move_id, still_ids
        self.move_name = f'part{move_id}'
        self.still_names = [f'part{still_id}' for still_id in still_ids]

        # collision check
        self.vertices_move = self.sim.get_body_vertices(self.move_name).T
        self.vertices_still = np.vstack([self.sim.get_body_vertices(still_name, world_frame=True).T for still_name in self.still_names])
        self.hull_move = trimesh.convex.convex_hull(self.vertices_move)
        self.hull_still = trimesh.convex.convex_hull(self.vertices_still)
        self.collision_manager = trimesh.collision.CollisionManager()
        self.collision_manager.add_object('hull_still', self.hull_still)
        
        # com
        coms = load_translation(assembly_dir)
        self.com_move = coms[move_id]
        self.coms_still = [coms[still_id] for still_id in self.still_ids]

        self.E0i_move = self.sim.get_body_E0i(self.move_name)
        self.E0is_still = [self.sim.get_body_E0i(still_name) for still_name in self.still_names]

        # state bounds
        self.min_box_move = self.vertices_move.min(axis=0)
        self.max_box_move = self.vertices_move.max(axis=0)
        self.size_box_move = self.max_box_move - self.min_box_move
        self.min_box_still = self.vertices_still.min(axis=0)
        self.max_box_still = self.vertices_still.max(axis=0)
        self.size_box_still = self.max_box_still - self.min_box_still
        self.state_lower_bound = (self.min_box_still - self.max_box_move) - 0.5 * self.size_box_move
        self.state_upper_bound = (self.max_box_still - self.min_box_move) + 0.5 * self.size_box_move

    def get_state(self):
        q = self.sim.get_joint_q(self.move_name)
        qdot = self.sim.get_joint_qdot(self.move_name)
        state = State(q, qdot)
        return state

    def set_state(self, state):
        # q, qdot = state.q, state.qdot
        q, qdot = state.q, np.zeros(self.ndof)
        self.sim.set_joint_state(self.move_name, q, qdot)

    def select_state(self, tree):
        raise NotImplementedError

    def select_action(self, tree, state):
        raise NotImplementedError

    def random_action(self):
        return np.random.random(self.ndof) * 2.0 - 1.0

    def apply_action(self, action):
        action = unit_vector(action) * self.force_mag # rotation, translation
        if len(action) == 3:
            force = np.concatenate([np.zeros(3), action])
        elif len(action) == 6:
            force = np.concatenate([action[:3] * 3, action[3:]])
        else:
            raise Exception
        self.sim.set_body_external_force(self.move_name, force)

    def q_distance(self, q0, q1):
        if self.rotation:
            boxes0 = transform_pts_by_state(np.vstack([self.min_box_move, self.max_box_move]), q0, com=self.com_move)
            boxes1 = transform_pts_by_state(np.vstack([self.min_box_move, self.max_box_move]), q1, com=self.com_move)
            return np.linalg.norm(boxes0 - boxes1, axis=1).sum()
        else:
            return np.linalg.norm(q0 - q1)

    def state_distance(self, state0, state1):
        return self.q_distance(state0.q, state1.q)

    def is_disassembled(self):
        E0i = self.sim.get_body_E0i(self.move_name)
        hull_move = self.hull_move.copy()
        hull_move.apply_transform(E0i)
        has_collision = self.collision_manager.in_collision_single(hull_move)
        if not has_collision:
            min_box_move, max_box_move = hull_move.vertices.min(axis=0), hull_move.vertices.max(axis=0)
            move_contain_still = (min_box_move <= self.min_box_still).all() and (max_box_move >= self.max_box_still).all()
            still_contain_move = (self.min_box_still <= min_box_move).all() and (self.max_box_still >= max_box_move).all()
            return not (move_contain_still or still_contain_move) # check if one hull fully contains another
        else:
            return False

    def get_path(self, tree, state):
        states = tree.get_root_path(state)
        path = []
        for state in states:
            path.append(state.q)
        return path

    def seed(self, seed):
        np.random.seed(seed)

    def render(self, path, record_path=None):
        if path is not None:
            self.sim.set_state_his(path, [np.zeros(self.ndof) for _ in range(len(path))])
        if record_path is None:
            SimRenderer.replay(self.sim)
        else:
            SimRenderer.replay(self.sim, record=True, record_path=record_path)

    def save_path(self, path, save_dir, n_save_state):
        save_path(save_dir, path, com=self.com_move, n_frame=n_save_state)

    def plan(self, max_time, seed=1, return_path=False, render=False, record_path=None):

        self.seed(seed)

        self.sim.reset()

        tree = Tree()
        init_state = self.get_state()
        tree.add_node(init_state)

        if self.is_disassembled():
            status = 'Start with goal'
            return (status, 0., None) if return_path else (status, 0.)

        status = 'Failure'
        path = None
        t_start = time()
        step = 0

        while True:

            if step == 0:
                state = init_state
                action = self.random_action()
            else:
                state = self.select_state(tree)
                action = self.select_action(tree, state)

            self.set_state(state)
            self.apply_action(action)
            self.sim.update_robot()

            states_between = []
            for _ in range(self.frame_skip):
                self.sim.forward(1)
                state_between = self.get_state()
                states_between.append(state_between)

                t_plan = time() - t_start
                if t_plan > max_time:
                    status = 'Timeout'
                    break

            new_state = states_between.pop()
            tree.add_node(new_state)
            tree.add_edge(state, new_state, action, states_between)

            if self.is_disassembled():
                status = 'Success'
                path = self.get_path(tree, new_state)
                break

            if status == 'Timeout':
                break

            step += 1

        if render:
            # tree.draw()
            self.render(path, record_path=record_path)

        return (status, t_plan, path) if return_path else (status, t_plan)

    def get_contact_bodies(self, part_id):
        return self.sim.get_contact_bodies(f'part{part_id}')


class BK_RRT(PhysicsPlanner):
    '''
    Behavioral Kinodynamic Rapidly-Exploring Random Trees (BK-RRT) [Zickler and Veloso 2009]
    '''
    def select_state(self, tree):

        # randomly sample a state (q)
        ratio = np.random.random(3)
        if self.rotation:
            translation = ratio * self.state_lower_bound + (1.0 - ratio) * self.state_upper_bound
            rotation = Rotation.random().as_rotvec()
            sampled_q = np.concatenate([translation, rotation])
        else:
            sampled_q = ratio * self.state_lower_bound + (1.0 - ratio) * self.state_upper_bound

        # find the closest state in the tree
        closest_state, closest_dist = None, np.inf
        for state in tree.get_nodes():
            dist = self.q_distance(sampled_q, state.q)
            if dist < closest_dist:
                closest_dist = dist
                closest_state = state

        return closest_state

    def select_action(self, tree, state):
        return self.random_action()


class BFSPlanner(PhysicsPlanner):

    # NOTE: can be tuned
    trans_dist_th = 0.05
    quat_dist_th = 0.5

    def get_trans_dist(self, state, new_state):
        return np.linalg.norm(state[:3] - new_state[:3])

    def get_quat(self, state):
        return Quaternion(Rotation.from_euler('xyz', state[3:]).as_quat()[[3, 0, 1, 2]])

    def get_quat_dist(self, state, new_state):
        return Quaternion.distance(self.get_quat(new_state), self.get_quat(state)) # NOTE: may use other distance functions

    def state_similar_r3(self, state, new_state):
        trans_dist = self.get_trans_dist(state, new_state)
        return trans_dist < self.trans_dist_th
        
    def state_similar_se3(self, state, new_state):
        trans_dist = self.get_trans_dist(state, new_state)
        quat_dist = self.get_quat_dist(state, new_state)
        return trans_dist < self.trans_dist_th and quat_dist < self.quat_dist_th

    def state_similar(self, state, new_state):
        if self.rotation:
            return self.state_similar_se3(state, new_state)
        else:
            return self.state_similar_r3(state, new_state)

    def any_state_similar(self, path, new_state):
        for state in path:
            if self.state_similar(state, new_state):
                return True
        return False

    def min_dist(self, path, new_state):
        if self.rotation:
            min_dist_trans, min_dist_rot = np.inf, np.inf
            for state in path:
                min_dist_trans = min(min_dist_trans, self.get_trans_dist(state, new_state))
                min_dist_rot = min(min_dist_rot, self.get_quat_dist(state, new_state))
            return min_dist_trans, min_dist_rot
        else:
            min_dist = np.inf
            for state in path:
                min_dist = min(min_dist, self.get_trans_dist(state, new_state))
            return min_dist

    def min_dist_separate(self, path, new_state):
        if self.rotation:
            min_dist_trans, cor_dist_rot = np.inf, np.inf
            cor_dist_trans, min_dist_rot = np.inf, np.inf
            for state in path:
                dist_trans = self.get_trans_dist(state, new_state)
                dist_rot = self.get_quat_dist(state, new_state)
                if dist_trans < min_dist_trans:
                    min_dist_trans = dist_trans
                    cor_dist_rot = dist_rot
                if dist_rot < min_dist_rot:
                    min_dist_rot = dist_rot
                    cor_dist_trans = dist_trans
            return (min_dist_trans, cor_dist_rot), (cor_dist_trans, min_dist_rot)
        else:
            min_dist = np.inf
            for state in path:
                min_dist = min(min_dist, self.get_trans_dist(state, new_state))
            return min_dist

    def random_rotate_actions(self, actions):
        random_vec = np.random.random(3)
        random_vec /= np.linalg.norm(random_vec)
        new_actions = []
        for action in actions:
            if actions.shape[1] == 3:
                new_action = np.cross(action, random_vec)
            elif actions.shape[1] == 6:
                new_action = np.concatenate([np.cross(action[:3], random_vec), np.cross(action[3:], random_vec)])
            else:
                raise Exception
            new_action /= np.linalg.norm(new_action)
            new_actions.append(new_action)
        return np.array(new_actions)

    def plan(self, *args, **kwargs):
        if self.rotation:
            return self.plan_rot(*args, **kwargs)
        else:
            return self.plan_trans(*args, **kwargs)

    def plan_trans(self, max_time, max_depth=None, seed=1, return_path=False, render=False, record_path=None):

        self.seed(seed)

        actions = np.array([
            [0, 0, 1], 
            [0, 0, -1],
            [0, 1, 0],
            [0, -1, 0],
            [1, 0, 0],
            [-1, 0, 0],
        ])
        # actions = self.random_rotate_actions(actions)

        status = 'Failure'
        path = None

        t_start = time()

        self.sim.reset()
        states = [[self.get_state(), []]]

        n_stages = 0

        while True: # stages

            state, curr_path = states.pop(0)

            for action in actions:

                temp_path = curr_path.copy()

                self.sim.reset()
                self.set_state(state)
                self.apply_action(action)

                while True:

                    self.set_state(self.get_state())

                    for _ in range(self.frame_skip):
                        self.sim.forward(1, verbose=False)
                        new_state = self.get_state()
                        temp_path.append(new_state.q)

                        t_plan = time() - t_start
                        if t_plan > max_time:
                            status = 'Timeout'
                            break

                    if self.is_disassembled():
                        status = 'Success'
                        path = temp_path
                        break

                    if status == 'Timeout':
                        break

                    if self.any_state_similar(temp_path[:-self.frame_skip], new_state.q):
                        break # back and forth

                if status in ['Success', 'Timeout']:
                    break

                states.append([new_state, temp_path])

            if status in ['Success', 'Timeout']:
                break

            n_stages += 1
            if n_stages == max_depth:
                break

        if render:
            self.render(path, record_path=record_path)

        return (status, t_plan, path) if return_path else (status, t_plan)

    def plan_rot(self, max_time, max_depth=None, seed=1, return_path=False, render=False, record_path=None):

        self.seed(seed)

        actions = np.array([
            [0, 0, 0, 0, 0, 1], 
            [0, 0, 0, 0, 0, -1],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, -1, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, -1, 0, 0],
            [0, 0, 1, 0, 0, 0], 
            [0, 0, -1, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, -1, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0],
            [-1, 0, 0, 0, 0, 0],
        ])
        # actions = self.random_rotate_actions(actions)

        status = 'Failure'
        path = None

        t_start = time()

        self.sim.reset()
        states = [[self.get_state(), []]]

        n_stages = 0

        while True: # stages

            stage_states = []

            while len(states) > 0:

                state, curr_path = states.pop(0)

                for action in actions:

                    temp_path = curr_path.copy()

                    self.sim.reset()
                    self.set_state(state)
                    self.apply_action(action)

                    while True:

                        self.set_state(self.get_state())

                        for _ in range(self.frame_skip):
                            self.sim.forward(1, verbose=False)
                            new_state = self.get_state()
                            temp_path.append(new_state.q)

                            t_plan = time() - t_start
                            if t_plan > max_time:
                                status = 'Timeout'
                                break

                        if self.is_disassembled():
                            status = 'Success'
                            path = temp_path
                            break

                        if status == 'Timeout':
                            break

                        if self.any_state_similar(temp_path[:-self.frame_skip], new_state.q):
                            # print(action, 'back and forth break')
                            break # back and forth

                    if status in ['Success', 'Timeout']:
                        break

                    stage_states.append([new_state, temp_path])

                    # if render:
                    #     print(f'Stage: {n_stages}, Action: {action}, Queue size: {len(states)}')
                    #     self.render(temp_path, record_path=record_path)

                if status in ['Success', 'Timeout']:
                    break

            if status in ['Success', 'Timeout']:
                break

            states = sorted(stage_states, key=lambda x: -len(x[1])) # sort based on path length
            n_stages += 1
            if n_stages == max_depth:
                break

        if render:
            self.render(path, record_path=record_path)

        return (status, t_plan, path) if return_path else (status, t_plan)


def get_planner(name):
    planners = {
        'bk-rrt': BK_RRT,
        'bfs': BFSPlanner,
    }
    assert name in planners, f'invalid planner name {name}'
    return planners[name]


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--planner', type=str, required=True, choices=['bfs', 'bk-rrt'])
    parser.add_argument('--id', type=str, required=True, help='assembly id (e.g. 00000)')
    parser.add_argument('--dir', type=str, default='joint_assembly', help='directory storing all assemblies')
    parser.add_argument('--move-id', type=int, default=0)
    parser.add_argument('--still-ids', type=int, nargs='+', default=[1])
    parser.add_argument('--rotation', default=False, action='store_true')
    parser.add_argument('--body-type', type=str, default='sdf', choices=['bvh', 'sdf'], help='simulation type of body')
    parser.add_argument('--sdf-dx', type=float, default=0.05, help='grid resolution of SDF')
    parser.add_argument('--collision-th', type=float, default=1e-2)
    parser.add_argument('--force-mag', type=float, default=100, help='magnitude of force')
    parser.add_argument('--frame-skip', type=int, default=100, help='control frequency')
    parser.add_argument('--max-time', type=float, default=120, help='timeout')
    parser.add_argument('--seed', type=int, default=1, help='random seed')
    parser.add_argument('--render', default=False, action='store_true', help='if render the result')
    parser.add_argument('--record-dir', type=str, default=None, help='directory to store rendering results')
    parser.add_argument('--save-dir', type=str, default=None)
    parser.add_argument('--n-save-state', type=int, default=100)
    parser.add_argument('--save-sdf', default=False, action='store_true')
    args = parser.parse_args()

    asset_folder = os.path.join(project_base_dir, './assets')
    assembly_dir = os.path.join(asset_folder, args.dir, args.id)

    if args.record_dir is None:
        record_path = None
    elif args.render:
        os.makedirs(args.record_dir, exist_ok=True)
        record_path = os.path.join(args.record_dir, args.id + '.gif')

    clear_saved_sdfs(assembly_dir)
    planner = get_planner(args.planner)(
        asset_folder, assembly_dir, args.move_id, args.still_ids, 
        args.rotation, args.body_type, args.sdf_dx, args.collision_th, args.force_mag, args.frame_skip, args.save_sdf
    )
    status, t_plan, path = planner.plan(
        args.max_time, seed=args.seed, return_path=True, render=args.render, record_path=record_path
    )
    clear_saved_sdfs(assembly_dir)

    print(f'Status: {status}, planning time: {t_plan}')

    if args.save_dir is not None:
        planner.save_path(path, args.save_dir, args.n_save_state)