import os
import sys

project_base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.append(project_base_dir)

import numpy as np
import redmax_py as redmax
import json
from argparse import ArgumentParser
from tqdm import tqdm

from assets.load import load_translation, load_part_ids
from assets.color import get_color
from utils.renderer import SimRenderer


def arr_to_str(arr):
    return ' '.join([str(x) for x in arr])


def get_xml_string(assembly_dir, force_part_id, body_type, sdf_dx, friction):
    part_ids = load_part_ids(assembly_dir)
    translation = load_translation(assembly_dir)
    if translation is None: translation = {part_id: [0, 0, 0] for part_id in part_ids}
    color_map = get_color(part_ids)
    body_type = body_type.upper()
    string = f'''
<redmax model="assemble">
<option integrator="BDF1" timestep="1e-3" gravity="0. 0. 1e-12"/>

<ground pos="0 0 -5" normal="0 0 1"/>
<default>
    <general_{body_type}_contact kn="1e6" kt="1e3" mu="{friction}" damping="0"/>
</default>
'''
    for part_id in ['0', '1']:
        joint_type = 'free3d-exp' if part_id == force_part_id else 'fixed'
        string += f'''
<robot>
    <link name="part{part_id}">
        <joint name="part{part_id}" type="{joint_type}" axis="0. 0. 0." pos="{arr_to_str(translation[part_id])}" quat="1 0 0 0" frame="WORLD" damping="0"/>
        <body name="part{part_id}" type="{body_type}" filename="{assembly_dir}/{part_id}.obj" pos="0 0 0" quat="1 0 0 0" scale="1 1 1" transform_type="OBJ_TO_JOINT" density="1" dx="{sdf_dx}" mu="0" rgba="{arr_to_str(color_map[part_id])}"/>
    </link>
</robot>
'''
    string += f'''
<contact>
'''
    string += f'''
    <general_{body_type}_contact general_body="part0" {body_type}_body="part1"/>
    <general_{body_type}_contact general_body="part1" {body_type}_body="part0"/>
'''
    string += f'''
</contact>
</redmax>
'''
    return string


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
    parser.add_argument('--dir', type=str, default='joint_assembly', help='directory storing all assemblies')
    parser.add_argument('--force', type=float, nargs='+', default=[0, 0, 0], help='external force to apply')
    parser.add_argument('--force-part-id', type=str, default='0', help='which part to apply force')
    parser.add_argument('--steps', type=int, default=1000, help='number of simulation steps')
    parser.add_argument('--body-type', type=str, default='sdf', choices=['bvh', 'sdf'], help='simulation type of body')
    parser.add_argument('--sdf-dx', type=float, default=0.05, help='grid resolution of SDF')
    parser.add_argument('--friction', type=float, default=0.0, help='friction coefficient between parts')
    parser.add_argument('--frame-skip', type=int, default=None)
    args = parser.parse_args()

    asset_folder = os.path.join(project_base_dir, './assets')
    assembly_dir = os.path.join(asset_folder, args.dir, args.id)

    assert args.force_part_id in ['0', '1']

    model_string = get_xml_string(
        assembly_dir=assembly_dir,
        force_part_id=args.force_part_id,
        body_type=args.body_type,
        sdf_dx=args.sdf_dx,
        friction=args.friction,
    )

    sim = redmax.Simulation(model_string, asset_folder)
    sim.reset()
    set_force(sim, args.force, args.force_part_id)

    for i in tqdm(range(args.steps)):
        sim.forward(1, verbose=True)

        if args.frame_skip is not None:
            if (i + 1) % args.frame_skip == 0:
                sim.zero_joint_qdot('part0')
        
    SimRenderer.replay(sim, record=False)
