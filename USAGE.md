# EDKMC Usage Guide

## 1. Installation

Install the calculator backend you plan to use before running the workflow.

### 1.1 MatPL

Official installation guide:

- `http://doc.lonxun.com/MatPL/install/Installation-online/`

After installation, check that one of these imports works:

```python
from src.ase.calculate import MatPL
```

or

```python
from src.ase.calculate import MatPL_calculator
```

This codebase will try `MatPL` first and then fall back to `MatPL_calculator`.

### 1.2 calorine

Official calculator page:

- `https://calorine.materialsmodeling.org/get_started/calculators.html`

Typical CPU-side installation:

```bash
pip install calorine --user
```

The calorine docs also show the CPU calculator as:

```python
from calorine.calculators import CPUNEP
calc = CPUNEP("nep.txt")
```

and the GPU calculator as:

```python
from calorine.calculators import GPUNEP
calc = GPUNEP("nep.txt")
```

### 1.3 Which backend to install

- Use MatPL if `calculator = matpl`
- Use calorine if `calculator = cpunep`
- If you may switch between them, install both

## 2. Files You Need

Required files:

- `edkmc.in`
- `edkmc.py` or `edkmc_pbc.py`
- the structure file given by `input_structure`
- the model file given by `model_file`

The input structure must contain:

- atomic coordinates
- a `group` array

The `group` array is used to:

- identify carbon atoms
- identify constrained atoms
- optionally identify metal atoms if `metal_elements` is not set

## 3. `edkmc.in`

Example:

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

Key parameters:

- `calculator`: backend selector, `matpl` or `cpunep`
- `model_file`: shared potential/model file
- `carbon_group`: group id of carbon atoms
- `constrain_group`: group id of fixed atoms during optimization
- `metal_elements`: preferred metal selector by element symbol
- `metal_group`: fallback metal selector if `metal_elements` is empty
- `cc_cutoff`: C-C neighbor cutoff
- `c_metal_cutoff`: C-metal neighbor cutoff
- `coord_threshold`: minimum C-C coordination number for an active carbon
- `fmax`: force threshold for BFGS
- `rotate_max`: maximum number of trial rotations
- `temperature`: temperature for Metropolis acceptance
- `early_stop`: stop after the first accepted move if `True`

## 4. Core Logic

Both scripts share the same workflow.

### `read_edkmc_input()`

Reads `edkmc.in`, converts types, and prints a parameter summary.

### `build_calculator()`

Builds the ASE calculator from `calculator` and `model_file`.

Supported backends:

- MatPL
- CPUNEP

### `optimize_structure()`

Runs constrained geometry optimization with ASE `BFGS`.

It:

- copies the structure
- attaches the calculator
- freezes atoms with `group == constrain_group`
- optimizes until `fmax` is reached

This function is used:

- once for the initial structure
- once after every trial rotation

### `find_active_carbon_atoms()`

Finds carbon atoms that satisfy both conditions:

- C-C coordination number is at least `coord_threshold`
- at least one metal atom is closer than `c_metal_cutoff`

Metal atoms are determined by:

- `metal_elements` first
- `metal_group` only as fallback

The code stops searching as soon as one valid metal neighbor is found.

### `rotate_carbon_bond()`

Creates one trial local bond-rotation move:

1. choose one active carbon `A`
2. find carbon neighbors of `A`
3. choose two neighbors `B` and `D`
4. define the local `A-B-D` plane
5. rotate the `A-B` pair by 90 degrees around the plane normal

### `monte_carlo_rotation()`

Drives the whole Monte Carlo loop:

1. find active carbon atoms
2. generate one trial rotation
3. optimize the rotated structure
4. compute `dE = E_new - E_old`
5. accept if `dE < 0`
6. otherwise accept with Metropolis probability

Accepted structures overwrite `model.xyz`.

Outputs:

- `model.xyz`
- `mc_energy.log`
- `mc_summary.log`
- `initial_optimization.traj`
- `initial_optimization.log`

## 5. Difference Between the Two Scripts

### `edkmc.py`

Use this for non-periodic or generic models.

### `edkmc_pbc.py`

Use this for periodic models.

The PBC version does two extra things:

- maps neighbor atoms to their nearest periodic images before rotation
- wraps the rotated structure back into the cell

That is the main difference between the two scripts.

## 6. How To Run

### Direct run

Non-periodic:

```bash
python edkmc.py
```

Periodic:

```bash
python edkmc_pbc.py
```

### Slurm run

The current `submit_job.sh` runs:

```bash
python edkmc.py
```

If your system is periodic, change that line to:

```bash
python edkmc_pbc.py
```

Then submit with:

```bash
sbatch submit_job.sh
```

## 7. `submit_job.sh` Workflow

The submission script does three things:

1. prepares the environment
2. runs repeated GPUMD growth/add-carbon steps
3. runs the Monte Carlo relaxation script

Current flow:

```text
10 x (GPUMD + add carbon)
-> 1 x edkmc.py
-> repeat 10 times
```

During each outer loop:

- `mc_energy.log` is archived to `mc_energyLog/`
- `mc_summary.log` is archived to `mc_summaryLog/`

## 8. Practical Run Order

1. Prepare `edkmc.in`
2. Make sure the structure has a `group` array
3. Decide whether the model is periodic
4. Use `edkmc.py` or `edkmc_pbc.py` accordingly
5. Test with a direct `python` run
6. Submit with `sbatch submit_job.sh`

## 9. Output Files

- `model.xyz`: latest accepted structure
- `mc_energy.log`: energy history during MC
- `mc_summary.log`: final statistics
- `initial_optimization.traj`: initial optimization trajectory
- `initial_optimization.log`: initial optimization log

## 10. Short Summary

- `edkmc.py` and `edkmc_pbc.py` share the same logic
- `edkmc_pbc.py` only changes the PBC rotation geometry
- `metal_elements` is the preferred way to define metals
- `submit_job.sh` should call `edkmc_pbc.py` for periodic systems
