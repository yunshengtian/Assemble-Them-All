# Assemble-Them-All

This repository contains the official code and dataset of [Assemble Them All: Physics-Based Planning for Generalizable Assembly by Disassembly (SIGGRAPH Asia 2022)](http://assembly.csail.mit.edu/). 

![representative](images/representative.gif)

**Authors**: Yunsheng Tian, Jie Xu, Yichen Li, Jieliang Luo, Shinjiro Sueda, Hui Li, Karl D.D. Willis, Wojciech Matusik

**Summary**: In this work, we propose a novel method to efficiently plan physically plausible assembly motion and sequences for real-world assemblies. Our method leverages the assembly-by-disassembly principle and physics-based simulation to efficiently explore a reduced search space. We also define a large-scale dataset consisting of thousands of physically valid industrial assemblies with a variety of assembly motions required.

## Installation

### 1. Clone repository (with dependencies)

```
git clone --recurse-submodules git@github.com:yunshengtian/Assemble-Them-All.git
```

### 2. Setup Python environment

```
conda env create -f environment.yml
conda activate assembly
```

**Note:** The package versions are critical to making sure the code runs smoothly.

### 3. Install simulation bindings

```
cd simulation
python setup.py install
```

**Note:** [CMake](https://cmake.org/) version > 3.16.0 is required.

### 4. Install assembly datasets

![dataset](images/dataset.png)

Our datasets are available at the links below (in .zip). Unzip them and place the folders under ``assets/`` of this code repository (recommended location, not required). 

[[Two-part main dataset (225.7MB)]](https://people.csail.mit.edu/yunsheng/Assemble-Them-All/joint_assembly.zip) [[Two-part rotational dataset (5.4MB)]](https://people.csail.mit.edu/yunsheng/Assemble-Them-All/joint_assembly_rotation.zip) [[Multi-part dataset (972.0MB)]](https://people.csail.mit.edu/yunsheng/Assemble-Them-All/multi_assembly.zip)

**Note 1:** Since a large portion of our dataset is modified from the [Fusion 360 Gallary Dataset](https://github.com/AutodeskAILab/Fusion360GalleryDataset), please refer to the [Fusion 360 Gallery Dataset License](https://github.com/AutodeskAILab/Fusion360GalleryDataset/blob/master/LICENSE.md) for legal usage.

**Note 2:** For point-based SDF collision check to work more accurately, we highly recommend subdividing the assembly meshes to have denser contact points by running ``assets/subdivide_batch.py`` (as we did in our experiments). For example, to subdivide the two-part assembly dataset saved in ``assets/joint_assembly`` and export to ``assets/joint_assembly_dense``:

```
python assets/subdivide_batch.py --source-dir assets/joint_assembly --target-dir assets/joint_assembly_dense --num-proc NUM_PROCESSES
```

To parallelize the subdivision, specify ``--num-proc`` with your preferred number of processes accordingly. After subdivision, please specify the new dataset directories for ``--dir`` argument to run the experiment scripts (see "Experiments" section below).

### 5. Test simulation

To verify if the installation steps are successful, run the following commands to test the two-part and multi-part assembly simulation respectively:

```
python examples/test_joint_sim.py --dir joint_assembly --id 00016 --gravity 9.8 --steps 2000
python examples/test_multi_sim.py --dir multi_assembly --id 00031 --gravity 9.8 --steps 2000
```

Then the simulation viewer will show up and produce the following animation:

|          00016 (two-part)           |         00031 (multi-part)          |
| :---------------------------------: | :---------------------------------: |
| ![00007](images/test_joint_sim.gif) | ![00007](images/test_multi_sim.gif) |

## Simulation Viewer

Some tips to interact with the simulation viewer:

- Scroll your mouse wheel to zoom in/out.
- Move your mouse while holding down the left mouse button to rotate the camera.
- Move your mouse while holding down ``Shift`` + left mouse button to move around the scene.
- Press ``S`` to slow down and ``F`` to speed up the simulation.
- Press ``Space`` to pause/resume the simulation.

## Experiments

### 1. Two-part assembly planning - single assembly

See below for commands and arguments for path planning on a single two-part assembly.

**Physics-based planning**

```
python examples/run_joint_plan.py --planner PLANNER --dir ASSEMBLY_DIR --id ASSEMBLY_ID
```

- ``--planner``: Path planning method. Choices: ``bfs`` (our method), ``bk-rrt``.

**Geometric-sampling-based planning**

```
python baselines/run_joint_plan.py --planner PLANNER --dir ASSEMBLY_DIR --id ASSEMBLY_ID
```

- ``--planner``: Path planning method. Choices: ``rrt``, ``rrt-connect``, ``birrt``, ``trrt``, ``matevec-trrt``.

**Other important arguments**

- ``--dir``: Directory of the two-part assembly dataset (e.g. ``joint_assembly`` if the dataset location is ``assets/joint_assembly``).
- ``--id``: Assembly id in the dataset (e.g. ``00000``).
- ``--rotation``: Add this argument if the assembly requires rotational motion.
- ``--max-time``: Planning timeout (in seconds, 120 by default).
- ``--render``: Add this argument to render the planned result.

For details on more arguments, please refer to the code. Below are some example animations that can be produced by the commands above.

|                   00007                   |                   00009                   |                   00016                   |
| :---------------------------------------: | :---------------------------------------: | :---------------------------------------: |
| ![00007](images/run_joint_plan/00007.gif) | ![00007](images/run_joint_plan/00009.gif) | ![00007](images/run_joint_plan/00016.gif) |

### 2. Two-part assembly planning - whole dataset

See below for commands and arguments for path planning on the whole two-part assembly dataset.

**Physics-based planning**

```
python examples/run_joint_plan_batch.py --planner PLANNER --dir ASSEMBLY_DIR --log-dir LOG_DIR
```

**Geometric-sampling-based planning**

```
python baselines/run_joint_plan_batch.py --planner PLANNER --dir ASSEMBLY_DIR --log-dir LOG_DIR
```

**Important arguments**

- ``--log-dir``: Directory to store the log files containing statistics of experiments.

For other arguments, see details in section "Two-part assembly planning - single assembly" above.

### 3. Multi-part assembly planning - single assembly

See below for commands and arguments for sequence (and path) planning on a single multi-part assembly.

**Physics-based planning**

```
python examples/run_multi_plan.py --seq-planner SEQ_PLANNER --path-planner PATH_PLANNER --dir ASSEMBLY_DIR --id ASSEMBLY_ID
```

- ``--seq-planner``: Sequence planning method. Choices: ``random``, ``queue``, ``prog-queue`` (our method).
- ``--path-planner``: Path planning method. Choices: ``bfs`` (our method), ``bk-rrt``.

**Geometric-sampling-based planning**

```
python baselines/run_multi_plan.py --seq-planner SEQ_PLANNER --path-planner PATH_PLANNER --dir ASSEMBLY_DIR --id ASSEMBLY_ID
```

- ``--seq-planner``: Sequence planning method. Choices: ``random``, ``queue`` (recommended).
- ``--path-planner``: Path planning method. Choices: ``rrt``, ``rrt-connect``, ``birrt``, ``trrt``, ``matevec-trrt``.

**Other important arguments**

- ``--dir``: Directory of the multi-part assembly dataset (e.g. ``multi_assembly`` if the dataset location is ``assets/multi_assembly``).
- ``--id``: Assembly id in the dataset (e.g. ``00000``).
- ``--rotation``: Add this argument if the assembly requires rotational motion.
- ``--seq-max-time``: Timeout for the whole sequence planning (in seconds, 3600 by default) with path planning in the loop.
- ``--path-max-time``: Timeout for the each path planning (in seconds, 120 by default).
- ``--render``: Add this argument to render the planned result.

For details on more arguments, please refer to the code. Below are some example animations that can be produced by the commands above.

|               00031 (step 1)                |               00031 (step 2)                |
| :-----------------------------------------: | :-----------------------------------------: |
| ![00007](images/run_multi_plan/00031_0.gif) | ![00007](images/run_multi_plan/00031_1.gif) |

### 4. Multi-part assembly planning - whole dataset

See below for commands and arguments for sequence (and path) planning on the whole multi-part assembly dataset.

**Physics-based planning**

```
python examples/run_multi_plan_batch.py --seq-planner SEQ_PLANNER --path-planner PATH_PLANNER --dir ASSEMBLY_DIR
```

**Geometric-sampling-based planning**

```
python baselines/run_multi_plan_batch.py --seq-planner SEQ_PLANNER --path-planner PATH_PLANNER --dir ASSEMBLY_DIR
```

**Important arguments**

- ``--log-dir``: Directory to store the log files containing statistics of experiments.

For other arguments, see details in section "Multi-part assembly planning - single assembly" above.

## Contact

Please feel free to contact yunsheng@csail.mit.edu or create a GitHub issue for any questions about the code or dataset.

## Citation

If you find our paper, code or dataset is useful, please consider citing:

```
@article{tian2022assemble,
    title={Assemble Them All: Physics-Based Planning for Generalizable Assembly by Disassembly},
    author={Tian, Yunsheng and Xu, Jie and Li, Yichen and Luo, Jieliang and Sueda, Shinjiro and Li, Hui and Willis, Karl D.D. and Matusik, Wojciech},
    journal={ACM Trans. Graph.},
    volume={41},
    number={6},
    articleno={278},
    numpages={15},
    year={2022},
    publisher={ACM}
}
```
