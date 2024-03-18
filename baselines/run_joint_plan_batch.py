import os
os.environ['OMP_NUM_THREADS'] = '1'
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

from argparse import ArgumentParser
import json

from utils.parallel import parallel_execute
from run_joint_plan import PyPlanner


def plan(assembly_dir, move_id, still_ids, rotation, planner_name, max_collision, adaptive_collision, step_size, max_time, seed, body_type, sdf_dx):
    planner = PyPlanner(assembly_dir, move_id, still_ids, rotation=rotation, max_collision=max_collision, adaptive_collision=adaptive_collision,
        body_type=body_type, sdf_dx=sdf_dx)
    status, t_plan = planner.plan(
        planner_name=planner_name, 
        step_size=step_size, 
        max_time=max_time, 
        seed=seed, 
        return_path=False,
        simplify=False, 
        render=False
    )
    id = os.path.basename(assembly_dir)
    return id, status, t_plan


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
    parser.add_argument('--dir', type=str, default='joint_assembly')
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

    worker_args = []
    for assembly_id in assembly_ids:
        assembly_dir = os.path.join(assemblies_dir, assembly_id)
        worker_args.append([assembly_dir, args.move_id, args.still_ids, args.rotation, args.planner, args.max_collision, args.adaptive_collision, args.step_size, args.max_time, args.seed, args.body_type, args.sdf_dx])

    result_status = {}
    count = 0

    try:
        for id, status, t_plan in parallel_execute(plan, worker_args, args.num_proc, show_progress=True):
            if status not in result_status:
                result_status[status] = {}
            result_status[status][id] = t_plan

            count += 1
            if count % 1 == 0:
                log_results(args.log_dir, result_status, verbose=False)

        log_results(args.log_dir, result_status, verbose=True)
    
    except KeyboardInterrupt:
        pass