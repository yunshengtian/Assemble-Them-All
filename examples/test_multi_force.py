import os
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

import numpy as np
import redmax_py as redmax
import json
from argparse import ArgumentParser
from tqdm import tqdm

from utils.renderer import SimRenderer
import time

def get_color(index):
    colors = [
        [210, 87, 89, 255],
        [237, 204, 73, 255],
        [60, 167, 221, 255],
        [190, 126, 208, 255],
        [108, 192, 90, 255],
    ]
    colors = np.array(colors) / 255.0
    return colors[int(index) % 5]


def arr_to_str(arr):
    return ' '.join([str(x) for x in arr])


def get_xml_string(assembly_dir, part_ids, force_part_id, body_type, sdf_dx, friction):
    with open(os.path.join(assembly_dir, 'translation.json'), 'r') as fp:
        translation = json.load(fp)
    body_type = body_type.upper()
    string = f'''
<redmax model="assemble">
<option integrator="BDF1" timestep="1e-3" gravity="0. 0. 1e-12"/>

<ground pos="0 0 -5" normal="0 0 1"/>
<default>
    <ground_contact kn="1e6" kt="1e3" mu="0.8" damping="5e1"/>
    <general_{body_type}_contact kn="1e6" kt="1e3" mu="{friction}" damping="0"/>
</default>
'''
    for part_id in part_ids:
        joint_type = 'free3d-exp' if part_id == force_part_id else 'fixed'
        string += f'''
<robot>
    <link name="part{part_id}">
        <joint name="part{part_id}" type="{joint_type}" axis="0. 0. 0." pos="{arr_to_str(translation[part_id])}" quat="1 0 0 0" frame="WORLD" damping="0"/>
        <body name="part{part_id}" type="{body_type}" filename="{assembly_dir}/{part_id}.obj" pos="0 0 0" quat="1 0 0 0" scale="1 1 1" transform_type="OBJ_TO_JOINT" density="1" dx="{sdf_dx}" mu="0" rgba="{arr_to_str(get_color(part_id))}"/>
    </link>
</robot>
'''
    string += f'''
<contact>
'''
    for part_id in part_ids:
        if part_id != force_part_id:
            string += f'''
    <general_{body_type}_contact general_body="part{part_id}" {body_type}_body="part{force_part_id}"/>
    <general_{body_type}_contact general_body="part{force_part_id}" {body_type}_body="part{part_id}"/>
'''
    string += f'''
</contact>
</redmax>
'''
    return string


def get_part_id(obj_name):
    return obj_name.split('_')[-1].replace('.obj', '')


def get_part_ids(obj_dir):
    part_ids = []
    for obj_name in os.listdir(obj_dir):
        if obj_name.endswith('.obj') and obj_name != 'assembly.obj':
            part_ids.append(get_part_id(obj_name))
    part_ids.sort(key=lambda x: int(x))
    return part_ids


def set_force(sim, force, part_id):
    '''
    Set force on one part
    '''
    if len(force) == 3: force = np.concatenate([np.zeros(3), np.array(force)])
    elif len(force) == 6: force = np.array(force)
    else: raise NotImplementedError
    if part_id is None: return
    sim.set_body_external_force(f'part{part_id}', force)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--id', type=str, required=True, help='assembly id (e.g. 00000)')
    parser.add_argument('--dir', type=str, default='multi_assembly', help='directory storing all assemblies')
    parser.add_argument('--part-ids', type=str, nargs='+', default=None, help='which parts to simulate (default: all)')
    parser.add_argument('--force', type=float, nargs='+', default=[0, 0, 0], help='external force to apply')
    parser.add_argument('--force-part-id', type=str, default=None, help='which part to apply force')
    parser.add_argument('--steps', type=int, default=1000, help='number of simulation steps')
    parser.add_argument('--body-type', type=str, default='sdf', choices=['bvh', 'sdf'], help='simulation type of body')
    parser.add_argument('--sdf-dx', type=float, default=0.05, help='grid resolution of SDF')
    parser.add_argument('--friction', type=float, default=0.0, help='friction coefficient between parts')
    parser.add_argument('--frame-skip', type=int, default=None)
    args = parser.parse_args()

    asset_dir = os.path.join(project_base_dir, './assets')
    assembly_dir = os.path.join(asset_dir, args.dir, args.id)

    if args.part_ids is None:
        part_ids = get_part_ids(assembly_dir)
    else:
        part_ids = args.part_ids
    if args.force_part_id is not None: assert args.force_part_id in part_ids

    model_string = get_xml_string(
        assembly_dir=assembly_dir,
        part_ids=part_ids,
        force_part_id=args.force_part_id,
        body_type=args.body_type,
        sdf_dx=args.sdf_dx,
        friction=args.friction,
    )

    sim = redmax.Simulation(model_string, asset_dir)
    sim.reset()
    set_force(sim, args.force, args.force_part_id)

    time_start = time.time()

    for i in tqdm(range(args.steps)):
        sim.forward(1, verbose=True)

        if args.frame_skip is not None:
            if (i + 1) % args.frame_skip == 0:
                sim.zero_joint_qdot(f'part{args.force_part_id}')
    
    time_end = time.time()

    print('time = ', time_end - time_start)

    SimRenderer.replay(sim, record=False)
