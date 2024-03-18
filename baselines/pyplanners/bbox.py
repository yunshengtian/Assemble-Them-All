import numpy as np
from scipy.spatial.transform import Rotation
from .transform import transform_pts_by_state


def transform_vectors(vectors, config):
    if len(config) == 3: # translation only
        return vectors + config
    elif len(config) == 6: # translation + rotation
        translation, rotation = config[:3], config[3:]
        rotation = Rotation.from_rotvec(rotation).as_matrix()
        return vectors.dot(rotation.T) + translation
    else:
        raise NotImplementedError


def get_bbox(vertices):
    return np.min(vertices, axis=0), np.max(vertices, axis=0)


def get_bbox_all(vertices, configs):
    min_box, max_box = np.zeros((0, 3)), np.zeros((0, 3))
    for config in configs:
        transformed_vertices = transform_pts_by_state(vertices, state=config)
        min_box = np.vstack([min_box, transformed_vertices.min(axis=0)])
        max_box = np.vstack([max_box, transformed_vertices.max(axis=0)])
    return min_box, max_box


def get_bbox_goal_test(vertices_move, vertices_still):
    min_box_move, max_box_move = get_bbox(vertices_move)
    min_box_still, max_box_still = get_bbox(vertices_still)
    
    def goal_test(config):
        if len(config) == 3: # translation only
            return ((config + min_box_move) >= max_box_still).any() or ((config + max_box_move) <= min_box_still).any()
        elif len(config) == 6: # translation + rotation
            vertices_move_new = transform_vectors(vertices_move, config)
            min_box_move_new, max_box_move_new = get_bbox(vertices_move_new)
            return (min_box_move_new >= max_box_still).any() or (max_box_move_new <= min_box_still).any()
        else:
            raise NotImplementedError
    return goal_test


def compute_nearest_bbox_goal(nodes, vertices_move, vertices_still):
    nodes_configs = np.array([n.config for n in nodes])
    if nodes_configs.shape[1] == 3: # translation only
        
        min_box_move, max_box_move = get_bbox(vertices_move)
        min_box_still, max_box_still = get_bbox(vertices_still)
        size_box_move = max_box_move - min_box_move
        dist_to_min = -(min_box_still - (nodes_configs + min_box_move) - size_box_move)
        dist_to_max = max_box_still - (nodes_configs + max_box_move) + size_box_move

    elif nodes_configs.shape[1] == 6: # translation + rotation

        min_box_move, max_box_move = get_bbox_all(vertices_move, nodes_configs)
        min_box_still, max_box_still = get_bbox(vertices_still)
        size_box_move = max_box_move - min_box_move
        dist_to_min = -(min_box_still - min_box_move - size_box_move)
        dist_to_max = max_box_still - max_box_move + size_box_move

    else:
        raise NotImplementedError
    
    assert (dist_to_min >= 0).all() and (dist_to_max >= 0).all()
    closest_to_min = np.unravel_index(np.argmin(dist_to_min), dist_to_min.shape)
    closest_to_max = np.unravel_index(np.argmin(dist_to_max), dist_to_max.shape)

    if dist_to_min[closest_to_min] < dist_to_max[closest_to_max]:
        node_idx, axis_idx = closest_to_min
        closest_node = nodes[node_idx]
        closest_goal = closest_node.config.copy()
        if nodes_configs.shape[1] == 3: # translation only
            closest_goal[axis_idx] = (min_box_still - min_box_move - size_box_move)[axis_idx]
        elif nodes_configs.shape[1] == 6: # translation + rotation
            closest_goal[axis_idx] += (min_box_still - min_box_move - size_box_move)[node_idx, axis_idx]
        else:
            raise NotImplementedError
    else:
        node_idx, axis_idx = closest_to_max
        closest_node = nodes[node_idx]
        closest_goal = closest_node.config.copy()
        if nodes_configs.shape[1] == 3: # translation only
            closest_goal[axis_idx] = (max_box_still - max_box_move + size_box_move)[axis_idx]
        elif nodes_configs.shape[1] == 6: # translation + rotation
            closest_goal[axis_idx] += (max_box_still - max_box_move + size_box_move)[node_idx, axis_idx]
        else:
            raise NotImplementedError

    return closest_goal, closest_node