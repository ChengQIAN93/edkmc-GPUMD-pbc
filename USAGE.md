# EDKMC Usage Guide

## 1. Installation

Before using `edkmc.py` or `edkmc_pbc.py`, install the calculator backend you plan to use.

### 1.1 Install MatPL

For MatPL, use the official online installation guide:

- MatPL installation guide: `http://doc.lonxun.com/MatPL/install/Installation-online/`

After installation, make sure one of the following imports works in your Python environment:

```python
from src.ase.calculate import MatPL
```

or

```python
from src.ase.calculate import MatPL_calculator
```

This codebase will first try `MatPL`, and then fall back to `MatPL_calculator` if needed.

### 1.2 Install calorine

According to the official calorine installation documentation, the stable version can be installed from PyPI with:

```bash
pip install calorine --user
```

The official documentation also notes that the PyPI package is distributed as source, so installation requires a C++11-compatible compiler such as GCC 4.8.1+ or Clang 3.3+. If the PyPI install fails, or if you need the latest development version, the official fallback command is:

```bash
pip install --user git+https://gitlab.com/materials-modeling/calorine.git
```

For calculator usage, the calorine documentation shows the CPU calculator as:

```python
from calorine.calculators import CPUNEP
calc = CPUNEP("nep.txt")
```

and the GPU calculator as:

```python
from calorine.calculators import GPUNEP
calc = GPUNEP("nep.txt")
```

### 1.3 Which backend should you install?

- If you want to run this workflow with `calculator = matpl`, install MatPL.
- If you want to run this workflow with `calculator = cpunep`, install calorine.
- If your environment may switch between the two, install both and choose the backend through `edkmc.in`.


## 2. Overview

This directory contains two Monte Carlo bond-rotation workflows:

- `edkmc.py`: standard workflow for non-periodic or general models
- `edkmc_pbc.py`: periodic version that applies nearest-image mapping before bond rotation

Both scripts use the same `edkmc.in` parameter file and the same overall logic:

1. Read parameters from `edkmc.in`
2. Read the input structure
3. Optimize the initial structure
4. Identify active carbon atoms
5. Perform a trial bond rotation
6. Re-optimize the rotated structure
7. Accept or reject the new structure using the Metropolis criterion
8. Repeat until `rotate_max` is reached, or stop early if `early_stop = True`


## 3. Required Files

Before running the workflow, make sure the following files are present:

- `edkmc.in`
- `edkmc.py` or `edkmc_pbc.py`
- structure file specified by `input_structure`
- potential/model file specified by `model_file`

The input structure must contain:

- atomic coordinates
- a `group` array

The `group` array is required because the code uses it to:

- identify the carbon group
- identify the constrained group
- optionally identify the metal group if `metal_elements` is not provided


## 4. Parameter File: `edkmc.in`

The parameter file is read by `read_edkmc_input()` in both scripts.

Current example:

```ini
calculator = matpl
model_file = nep13B.txt

carbon_group = 1
constrain_group = 0
metal_group = 2
metal_elements = Cu

cc_cutoff = 1.7
c_metal_cutoff = 4.5
coord_threshold = 2

fmax = 0.05

input_structure = restart.xyz
rotate_max = 30
temperature = 1400.0
early_stop = False

optimized_file = optimized.xyz
final_structure = model.xyz
```

### Key parameters

- `calculator`
  - `matpl`: use GPU MatPL
  - `cpunep`: use CPU CPUNEP

- `model_file`
  - shared model/potential file for both calculators

- `carbon_group`
  - atoms in this group are treated as carbon atoms for the MC logic

- `constrain_group`
  - atoms in this group are fixed during geometry optimization

- `metal_elements`
  - preferred way to identify metal atoms
  - examples:
    - `metal_elements = Cu`
    - `metal_elements = Cu Zn`
    - `metal_elements = Cu,Zn`

- `metal_group`
  - fallback metal definition if `metal_elements` is empty

- `cc_cutoff`
  - cutoff for carbon-carbon neighbors

- `c_metal_cutoff`
  - cutoff for carbon-metal neighbors

- `coord_threshold`
  - minimum C-C coordination number required for a carbon atom to be considered active

- `fmax`
  - force convergence threshold for BFGS optimization

- `rotate_max`
  - maximum number of Monte Carlo trial rotations

- `temperature`
  - temperature used in the Metropolis acceptance probability

- `early_stop`
  - if `True`, the MC loop stops after the first accepted move


## 5. Core Workflow Logic

### 5.1 Parameter parsing

Function:

- `read_edkmc_input()` in [edkmc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc.py:11)
- `read_edkmc_input()` in [edkmc_pbc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc_pbc.py:11)

What it does:

- reads `edkmc.in`
- converts integers, floats, booleans, and string/list-like parameters
- prints a parameter summary before the simulation starts


### 5.2 Calculator construction

Function:

- `build_calculator()` in [edkmc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc.py:103)
- `build_calculator()` in [edkmc_pbc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc_pbc.py:103)

What it does:

- builds the ASE calculator from `calculator` and `model_file`
- supports:
  - `MatPL`
  - `CPUNEP`

This is the place where the model backend is selected.


### 5.3 Structure optimization

Function:

- `optimize_structure()` in [edkmc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc.py:134)
- `optimize_structure()` in [edkmc_pbc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc_pbc.py:134)

What it does:

- copies the input structure
- attaches the calculator
- constrains all atoms with `group == constrain_group`
- runs ASE `BFGS`
- returns:
  - optimized structure
  - final energy

This function is used twice:

- once for the initial structure
- once after every trial bond rotation


### 5.4 Active carbon detection

Function:

- `find_active_carbon_atoms()` in [edkmc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc.py:182)
- `find_active_carbon_atoms()` in [edkmc_pbc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc_pbc.py:182)

This is one of the most important functions in the workflow.

It identifies a carbon atom as active only if both conditions are satisfied:

1. its C-C coordination number is at least `coord_threshold`
2. it has at least one metal neighbor within `c_metal_cutoff`

Important implementation details:

- carbon atoms are defined by `carbon_group`
- metal atoms are defined by:
  - `metal_elements` first
  - `metal_group` only as fallback
- the code only checks metal neighbors for carbons that already pass the C-C coordination threshold
- once one valid metal neighbor is found, the search stops immediately

This makes the logic both more physically controlled and more efficient.


### 5.5 Trial bond rotation

Function:

- `rotate_carbon_bond()` in [edkmc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc.py:294)
- `rotate_carbon_bond()` in [edkmc_pbc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc_pbc.py:300)

What it does:

1. randomly choose one active carbon atom `A`
2. find carbon neighbors around `A`
3. randomly choose two neighbors `B` and `D`
4. define a local plane `A-B-D`
5. rotate the `A-B` pair by 90 degrees around the plane normal

### Non-periodic version

`edkmc.py` rotates directly using the current Cartesian positions.

### Periodic version

`edkmc_pbc.py` first maps neighbors to their nearest periodic images before building the rotation geometry.

Related helper:

- `map_to_nearest_image()` in [edkmc_pbc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc_pbc.py:273)

That is the core difference between the two programs.


### 5.6 Monte Carlo acceptance

Function:

- `monte_carlo_rotation()` in [edkmc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc.py:406)
- `monte_carlo_rotation()` in [edkmc_pbc.py](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/edkmc_pbc.py:444)

This is the main driver of the simulation.

For each MC step:

1. find active carbon atoms
2. generate one trial rotation
3. optimize the rotated structure
4. compute `delta_E = E_new - E_old`
5. accept or reject:
   - always accept if `E_new < E_old`
   - otherwise accept with probability:

```text
P = exp(-delta_E / (kB * T))
```

If a move is accepted:

- the current structure is updated
- `model.xyz` is overwritten with the accepted structure

If a move is rejected:

- the previous accepted structure is kept

Outputs:

- `model.xyz`
- `mc_energy.log`
- `mc_summary.log`
- `initial_optimization.traj`
- `initial_optimization.log`


## 6. Which Script Should Be Used?

### Use `edkmc.py` if

- your structure is non-periodic
- or you do not need nearest-image handling during the bond rotation step

### Use `edkmc_pbc.py` if

- your structure is periodic
- bond neighbors may cross cell boundaries
- the local rotation geometry should be built using nearest periodic images

For periodic slab or bulk-like models, `edkmc_pbc.py` is usually the correct choice.


## 7. Direct Command-Line Usage

### Non-periodic workflow

```bash
python edkmc.py
```

### Periodic workflow

```bash
python edkmc_pbc.py
```

Both commands assume:

- the current working directory contains `edkmc.in`
- `input_structure` exists
- `model_file` exists


## 8. How `submit_job.sh` Works

The current submission script is [submit_job.sh](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/submit_job.sh:1).

Its logic is:

1. request 1 GPU job from Slurm
2. load compiler and CUDA modules
3. activate the `calorine` conda environment
4. create output directories:
   - `dump`
   - `thermo`
   - `movie`
   - `mc_energyLog`
   - `mc_summaryLog`
5. run two nested loops:
   - inner loop `j = 1..10`
   - outer loop `i = 1..10`

### Inner loop meaning

In each inner-loop iteration, the script does:

1. run GPUMD
2. remove old `model.xyz`
3. run `python add-carbon-gpumd.py`
4. archive:
   - `dump.xyz`
   - `thermo.out`
   - `movie.xyz`
5. delete the temporary GPUMD outputs

This part appears to generate growth trajectories and intermediate structures.

### Outer loop meaning

After 10 GPUMD/add-carbon cycles, the script runs:

```bash
python edkmc.py
```

Then it archives:

- `mc_energy.log`
- `mc_summary.log`

So the current job flow is:

```text
10 x (GPUMD + add carbon)
-> 1 x EDKMC optimization
-> repeat 10 times
```


## 9. Important Note About `submit_job.sh`

Right now, [submit_job.sh](/d:/BaiduSyncdisk/edkmc-gpumd/kmc/submit_job.sh:29) runs:

```bash
python edkmc.py
```

If your model is periodic and you want to use the periodic-safe rotation logic, you should change that line to:

```bash
python edkmc_pbc.py
```

That is the single most important usage change for periodic systems.


## 10. Recommended Usage Patterns

### Case A: non-periodic model

Use:

```bash
python edkmc.py
```

and keep in `submit_job.sh`:

```bash
python edkmc.py
```

### Case B: periodic model

Use:

```bash
python edkmc_pbc.py
```

and modify `submit_job.sh`:

```bash
python edkmc_pbc.py
```


## 11. Suggested Run Procedure

### Step 1

Prepare `edkmc.in`

Make sure these are correct:

- `calculator`
- `model_file`
- `input_structure`
- `carbon_group`
- `constrain_group`
- `metal_elements` or `metal_group`

### Step 2

Make sure the input structure has a valid `group` array.

### Step 3

Decide whether the model is periodic:

- non-periodic: use `edkmc.py`
- periodic: use `edkmc_pbc.py`

### Step 4

Run directly for testing:

```bash
python edkmc.py
```

or

```bash
python edkmc_pbc.py
```

### Step 5

After local validation, submit through Slurm:

```bash
sbatch submit_job.sh
```

### Step 6

Check outputs:

- Slurm stdout/stderr
- `mc_energy.log`
- `mc_summary.log`
- `model.xyz`


## 12. Practical Interpretation of Outputs

- `model.xyz`
  - the latest accepted structure

- `mc_energy.log`
  - step-by-step MC energy history
  - useful for checking whether moves are being accepted

- `mc_summary.log`
  - final summary of the MC run
  - useful for acceptance statistics

- `initial_optimization.traj`
  - the first optimization trajectory

- `initial_optimization.log`
  - the first optimization log


## 13. Minimal Example

For a periodic Cu-supported model, a typical setup is:

```ini
calculator = matpl
model_file = nep13B.txt

carbon_group = 1
constrain_group = 0
metal_group = 2
metal_elements = Cu

cc_cutoff = 1.7
c_metal_cutoff = 4.5
coord_threshold = 2

fmax = 0.05

input_structure = restart.xyz
rotate_max = 30
temperature = 1400.0
early_stop = False
```

Then:

```bash
python edkmc_pbc.py
```

If submitted through Slurm, update `submit_job.sh` so that it also calls:

```bash
python edkmc_pbc.py
```


## 14. Short Summary

If you only remember three things:

1. `edkmc.py` is the standard version, `edkmc_pbc.py` is the periodic version.
2. `metal_elements` is the preferred way to define metal atoms.
3. For periodic models, `submit_job.sh` should run `python edkmc_pbc.py`, not `python edkmc.py`.
