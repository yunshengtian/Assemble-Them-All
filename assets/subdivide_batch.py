'''
Batch subdivision for meshes
'''

import os
os.environ['OMP_NUM_THREADS'] = '1'
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

from argparse import ArgumentParser

from utils.parallel import parallel_execute
from subdivide import subdivide_assembly


parser = ArgumentParser()
parser.add_argument('--source-dir', type=str, required=True)
parser.add_argument('--target-dir', type=str, required=True)
parser.add_argument('--max-edge', type=float, default=0.5)
parser.add_argument('--num-proc', type=int, default=8)
args = parser.parse_args()

source_ids = []
for source_id in os.listdir(args.source_dir):
    dir_path = os.path.join(args.source_dir, source_id)
    if os.path.isdir(dir_path):
        source_ids.append(source_id)
source_ids.sort()

worker_args = []
for source_id in source_ids:
    worker_args.append([
        os.path.join(args.source_dir, source_id),
        os.path.join(args.target_dir, source_id),
        args.max_edge,
    ])

os.makedirs(args.target_dir, exist_ok=True)
for _ in parallel_execute(subdivide_assembly, worker_args, args.num_proc):
    pass
