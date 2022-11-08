import os
os.environ['OMP_NUM_THREADS'] = '1'
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

from argparse import ArgumentParser
import json

from assets.save import clear_saved_sdfs
from utils.parallel import parallel_execute
from run_multi_plan import get_seq_planner


def seq_plan(assembly_dir, seq_planner_name, path_planner_name, 
    rotation, max_collision, adaptive_collision, step_size, seq_max_time, path_max_time, seed, body_type, sdf_dx,
    save_dir, n_save_state):

    clear_saved_sdfs(assembly_dir)
    seq_planner = get_seq_planner(seq_planner_name)(assembly_dir)
    seq_status, sequence, seq_count, t_plan = seq_planner.plan_sequence(path_planner_name, 
        rotation, body_type, sdf_dx, max_collision, adaptive_collision, step_size, seq_max_time, path_max_time, seed, 
        render=False, save_dir=save_dir, n_save_state=n_save_state, verbose=False)
    clear_saved_sdfs(assembly_dir)

    id = os.path.basename(assembly_dir)
    return id, seq_status, sequence, seq_count, t_plan


def log_results(log_dir, result_status, verbose=False):
    stats = []
    for key in result_status.keys():
        key_stats = f'{key}: {len(result_status[key])}'
        if verbose:
            print(key_stats)
        stats.append(key_stats)
        with open(os.path.join(log_dir, f'{key}.json'), 'w') as fp:
            json.dump(dict(sorted(result_status[key].items())), fp)
    with open(os.path.join(log_dir, 'stats.txt'), 'w') as fp:
        fp.writelines([x + '\n' for x in stats])


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--dir', type=str, default='multi_assembly', help='directory storing all assemblies')
    parser.add_argument('--rotation', default=False, action='store_true')
    parser.add_argument('--seq-planner', type=str, required=True, choices=['random', 'queue'])
    parser.add_argument('--path-planner', type=str, required=True, choices=['rrt', 'rrt-connect', 'birrt', 'trrt', 'matevec-trrt'])
    parser.add_argument('--body-type', type=str, default='sdf', choices=['bvh', 'sdf'], help='simulation type of body')
    parser.add_argument('--sdf-dx', type=float, default=0.05, help='grid resolution of SDF')
    parser.add_argument('--max-collision', type=float, default=1e-2)
    parser.add_argument('--adaptive-collision', default=False, action='store_true')
    parser.add_argument('--step-size', type=float, default=1e-2)
    parser.add_argument('--seq-max-time', type=float, default=3600, help='sequence planning timeout')
    parser.add_argument('--path-max-time', type=float, default=120, help='path planning timeout')
    parser.add_argument('--seed', type=int, default=1, help='random seed')
    parser.add_argument('--render', default=False, action='store_true', help='if render the result')
    parser.add_argument('--save-dir', type=str, default=None)
    parser.add_argument('--n-save-state', type=int, default=100)
    parser.add_argument('--num-proc', type=int, default=8)
    parser.add_argument('--log-dir', type=str, required=True)
    args = parser.parse_args()

    asset_folder = os.path.join(project_base_dir, './assets')
    assemblies_dir = os.path.join(asset_folder, args.dir)
    assembly_ids = []
    for assembly_id in os.listdir(assemblies_dir):
        assembly_dir = os.path.join(assemblies_dir, assembly_id)
        if os.path.isdir(assembly_dir):
            assembly_ids.append(assembly_id)
    assembly_ids.sort()

    os.makedirs(args.log_dir)
    with open(os.path.join(args.log_dir, 'config.json'), 'w') as fp:
        json.dump(vars(args), fp)

    if args.rotation: args.seq_max_time *= 2

    worker_args = []
    for assembly_id in assembly_ids:
        assembly_dir = os.path.join(assemblies_dir, assembly_id)
        worker_args.append([
            assembly_dir, args.seq_planner, args.path_planner,
            args.rotation, args.max_collision, args.adaptive_collision, args.step_size, args.seq_max_time, args.path_max_time, args.seed, args.body_type, args.sdf_dx,
            args.save_dir, args.n_save_state
        ])

    result_status = {}
    count = 0

    try:
        for id, seq_status, sequence, seq_count, t_plan in parallel_execute(seq_plan, worker_args, args.num_proc, show_progress=True):
            if seq_status not in result_status:
                result_status[seq_status] = {}
            result_status[seq_status][id] = [seq_count, t_plan, sequence]

            count += 1
            if count % 1 == 0:
                log_results(args.log_dir, result_status, verbose=False)

        log_results(args.log_dir, result_status, verbose=True)

    except KeyboardInterrupt:
        pass