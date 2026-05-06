import random

import numpy as np
from ase.constraints import FixAtoms
from ase.geometry import get_distances
from ase.io import read, write
from ase.optimize import BFGS
from scipy.spatial.distance import cdist


def read_edkmc_input(filename="edkmc.in"):
    """
    Read simulation parameters from an edkmc.in file.

    Parameters
    ----------
    filename : str
        Input parameter file name.

    Returns
    -------
    dict
        Dictionary containing all parsed parameters.
    """
    params = {
        # Structure optimization parameters
        "constrain_group": 2,
        "fmax": 0.05,
        "calculator": "matpl",
        "model_file": "nep.txt",

        # Active carbon identification parameters
        "carbon_group": 1,
        "metal_group": 2,
        "metal_elements": [],
        "cc_cutoff": 1.7,
        "c_metal_cutoff": 3.0,
        "coord_threshold": 3,

        # Monte Carlo parameters
        "input_structure": "model.xyz",
        "rotate_max": 10,
        "temperature": 300.0,
        "early_stop": True,

        # Output files
        "optimized_file": "optimized.xyz",
        "final_structure": "model.xyz",
    }

    try:
        with open(filename, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                if "=" not in line:
                    continue

                key, value = line.split("=", 1)
                key = key.strip()
                value = value.split("#")[0].strip()

                if key in ["constrain_group", "carbon_group", "metal_group", "coord_threshold", "rotate_max"]:
                    params[key] = int(value)
                elif key in ["fmax", "cc_cutoff", "c_metal_cutoff", "temperature"]:
                    params[key] = float(value)
                elif key in ["early_stop"]:
                    params[key] = value.lower() in ["true", "yes", "1"]
                elif key == "metal_elements":
                    params[key] = [item.strip() for item in value.replace(",", " ").split() if item.strip()]
                elif key == "nep_file":
                    params["model_file"] = value
                elif key in [
                    "optimized_file",
                    "input_structure",
                    "final_structure",
                    "calculator",
                    "model_file",
                ]:
                    params[key] = value

        print(f"Successfully read parameter file: {filename}")
    except FileNotFoundError:
        print(f"Warning: parameter file {filename} was not found. Using default values.")

    print("\n" + "=" * 70)
    print("Monte Carlo Parameter Summary")
    print("=" * 70)
    print(f"Input structure:       {params['input_structure']}")
    print(f"Calculator:            {params['calculator']}")
    print(f"Model file:            {params['model_file']}")
    print(f"Max rotation steps:    {params['rotate_max']}")
    print(f"Temperature:           {params['temperature']} K")
    print(f"Early stop:            {params['early_stop']}")
    print(f"Metal elements:        {params['metal_elements']}")
    print("=" * 70)

    return params


def build_calculator(params):
    """
    Build an ASE calculator from the input parameters.
    """
    calculator_name = str(params.get("calculator", "matpl")).strip().lower()

    if calculator_name in ["matpl", "gpu", "gpu_matpl"]:
        try:
            from src.ase.calculate import MatPL as MatPLCalculator
        except ImportError:
            from src.ase.calculate import MatPL_calculator as MatPLCalculator

        model_file = params.get("model_file")
        if not model_file:
            raise ValueError("MATPL requires 'model_file'.")
        return MatPLCalculator(model_file=model_file)

    if calculator_name in ["cpunep", "cpu", "calorine"]:
        from calorine.calculators import CPUNEP

        model_file = params.get("model_file")
        if not model_file:
            raise ValueError("CPUNEP requires 'model_file'.")
        return CPUNEP(model_file)

    raise ValueError(
        f"Unsupported calculator '{params.get('calculator')}'. "
        "Use 'matpl' or 'cpunep'."
    )


def optimize_structure(atoms, params, silent=True):
    """
    Minimize the structure energy while constraining one atom group.

    Parameters
    ----------
    atoms : ASE Atoms
        Input structure.
    params : dict
        Simulation parameter dictionary.
    silent : bool
        If True, do not write optimizer log or trajectory files.

    Returns
    -------
    tuple
        Optimized atoms object and final energy in eV.
    """
    atoms_opt = atoms.copy()

    calc = build_calculator(params)
    atoms_opt.calc = calc

    if "group" not in atoms_opt.arrays:
        raise ValueError("The input structure does not contain a 'group' array.")

    groups = atoms_opt.arrays["group"]
    constrain_group = params["constrain_group"]
    group_indices = np.where(groups == constrain_group)[0]

    if len(group_indices) > 0:
        atoms_opt.set_constraint(FixAtoms(indices=group_indices))

    if silent:
        optimizer = BFGS(atoms_opt, logfile=None)
    else:
        optimizer = BFGS(
            atoms_opt,
            trajectory="initial_optimization.traj",
            logfile="initial_optimization.log",
        )

    optimizer.run(fmax=params["fmax"])
    final_energy = atoms_opt.get_potential_energy()

    return atoms_opt, final_energy


def find_active_carbon_atoms(atoms, params, verbose=True):
    """
    Identify active carbon atoms.

    A carbon atom is considered active when:
    1. Its C-C coordination number is at least coord_threshold.
    2. It has at least one nearby metal atom within c_metal_cutoff.

    Metal atoms are selected by element symbol if metal_elements is provided.
    Otherwise, the function falls back to metal_group.

    Parameters
    ----------
    atoms : ASE Atoms
        Input structure.
    params : dict
        Simulation parameter dictionary.
    verbose : bool
        If True, print summary information.

    Returns
    -------
    numpy.ndarray
        Indices of active carbon atoms.
    """
    groups = atoms.arrays["group"]
    positions = atoms.get_positions()

    carbon_group = params["carbon_group"]
    cc_cutoff = params["cc_cutoff"]
    c_metal_cutoff = params["c_metal_cutoff"]
    coord_threshold = params["coord_threshold"]

    raw_metal_elements = params.get("metal_elements", [])
    if isinstance(raw_metal_elements, str):
        metal_elements = {
            item.strip()
            for item in raw_metal_elements.replace(",", " ").split()
            if item.strip()
        }
    else:
        metal_elements = {str(item).strip() for item in raw_metal_elements if str(item).strip()}

    carbon_mask = groups == carbon_group
    carbon_indices = np.where(carbon_mask)[0]
    num_carbon = len(carbon_indices)
    symbols = np.array(atoms.get_chemical_symbols())

    if metal_elements:
        metal_mask = np.isin(symbols, list(metal_elements))
    else:
        metal_group = params["metal_group"]
        metal_mask = groups == metal_group

    metal_mask = metal_mask & (~carbon_mask)
    metal_indices = np.where(metal_mask)[0]
    use_pbc = np.any(atoms.get_pbc())

    # Count C-C neighbors for each carbon atom.
    if use_pbc:
        cc_distances = np.zeros((num_carbon, num_carbon))
        for i, idx_i in enumerate(carbon_indices):
            for j, idx_j in enumerate(carbon_indices):
                if i == j:
                    continue
                _, dist = get_distances(
                    [positions[idx_i]],
                    [positions[idx_j]],
                    cell=atoms.get_cell(),
                    pbc=atoms.get_pbc(),
                )
                cc_distances[i, j] = dist[0, 0]
    else:
        cc_distances = cdist(positions[carbon_indices], positions[carbon_indices])

    cc_neighbor_count = np.sum((cc_distances > 0) & (cc_distances < cc_cutoff), axis=1)

    # Only inspect metal neighbors for carbons that already pass the C-C threshold.
    # Stop searching as soon as one valid metal neighbor is found.
    has_metal_neighbor = np.zeros(num_carbon, dtype=bool)
    candidate_local_indices = np.where(cc_neighbor_count >= coord_threshold)[0]

    if len(metal_indices) > 0:
        for local_idx in candidate_local_indices:
            carbon_idx = carbon_indices[local_idx]
            carbon_position = positions[carbon_idx]

            for metal_idx in metal_indices:
                if use_pbc:
                    _, dist = get_distances(
                        [carbon_position],
                        [positions[metal_idx]],
                        cell=atoms.get_cell(),
                        pbc=atoms.get_pbc(),
                    )
                    distance = dist[0, 0]
                else:
                    distance = np.linalg.norm(carbon_position - positions[metal_idx])

                if distance < c_metal_cutoff:
                    has_metal_neighbor[local_idx] = True
                    break

    active_mask = (cc_neighbor_count >= coord_threshold) & has_metal_neighbor
    active_indices = carbon_indices[np.where(active_mask)[0]]

    if verbose:
        print(f"Found {len(active_indices)} active carbon atoms.")

    return active_indices


def rotate_carbon_bond(atoms, active_indices, params, verbose=True):
    """
    Perform a trial rotation on a local carbon-carbon bond environment.

    The procedure is:
    1. Randomly choose one active carbon atom A.
    2. Randomly choose two distinct carbon neighbors B and D around A.
    3. Rotate the A-B pair by 90 degrees around the normal of the ABD plane.

    Parameters
    ----------
    atoms : ASE Atoms
        Input structure.
    active_indices : array-like
        Indices of active carbon atoms.
    params : dict
        Simulation parameter dictionary.
    verbose : bool
        If True, print rotation details.

    Returns
    -------
    tuple
        Rotated atoms object and a success flag.
    """
    if len(active_indices) == 0:
        if verbose:
            print("Error: no active carbon atoms are available for selection.")
        return atoms, False

    atom_A_idx = np.random.choice(active_indices)

    carbon_group = params["carbon_group"]
    cc_cutoff = params["cc_cutoff"]
    groups = atoms.arrays["group"]
    positions = atoms.get_positions()

    carbon_indices = np.where(groups == carbon_group)[0]
    pos_A = positions[atom_A_idx]
    neighbor_candidates_A = []

    for idx in carbon_indices:
        if idx == atom_A_idx:
            continue

        if np.any(atoms.get_pbc()):
            _, dist = get_distances(
                [pos_A],
                [positions[idx]],
                cell=atoms.get_cell(),
                pbc=atoms.get_pbc(),
            )
            distance = dist[0, 0]
        else:
            distance = np.linalg.norm(positions[idx] - pos_A)

        if distance < cc_cutoff:
            neighbor_candidates_A.append((idx, distance))

    if len(neighbor_candidates_A) < 2:
        if verbose:
            print(f"Error: atom A (index {atom_A_idx}) has fewer than two carbon neighbors.")
        return atoms, False

    selected_indices = np.random.choice(len(neighbor_candidates_A), size=2, replace=False)
    atom_B_idx, dist_AB = neighbor_candidates_A[selected_indices[0]]
    atom_D_idx, dist_AD = neighbor_candidates_A[selected_indices[1]]

    if verbose:
        print(f"  Selected active carbon A: index {atom_A_idx}")
        print(f"  Neighbor carbon B: index {atom_B_idx}, distance {dist_AB:.4f} A")
        print(f"  Neighbor carbon D: index {atom_D_idx}, distance {dist_AD:.4f} A")

    pos_B = positions[atom_B_idx]
    pos_D = positions[atom_D_idx]

    vec_AB = pos_B - pos_A
    vec_AD = pos_D - pos_A
    normal_vector = np.cross(vec_AB, vec_AD)

    if np.linalg.norm(normal_vector) < 1e-6:
        if verbose:
            print("Error: atoms A, B, and D are collinear, so the rotation plane is undefined.")
        return atoms, False

    rotation_axis = normal_vector / np.linalg.norm(normal_vector)
    rotation_center = (pos_A + pos_B) / 2.0
    theta = np.pi / 2.0

    def rodrigues_rotation(point, axis, center, angle):
        p = point - center
        cos_theta = np.cos(angle)
        sin_theta = np.sin(angle)
        p_rot = (
            p * cos_theta
            + np.cross(axis, p) * sin_theta
            + axis * np.dot(axis, p) * (1 - cos_theta)
        )
        return p_rot + center

    new_pos_A = rodrigues_rotation(pos_A, rotation_axis, rotation_center, theta)
    new_pos_B = rodrigues_rotation(pos_B, rotation_axis, rotation_center, theta)

    atoms_rotated = atoms.copy()
    new_positions = atoms_rotated.get_positions()
    new_positions[atom_A_idx] = new_pos_A
    new_positions[atom_B_idx] = new_pos_B
    atoms_rotated.set_positions(new_positions)

    return atoms_rotated, True


def monte_carlo_rotation(params):
    """
    Run the Monte Carlo bond rotation workflow.

    Parameters
    ----------
    params : dict
        Simulation parameter dictionary.

    Returns
    -------
    dict
        Final simulation results.
    """
    kB = 8.617333262e-5

    print("\n" + "#" * 70)
    print("#" + " " * 68 + "#")
    print("#" + " " * 18 + "Monte Carlo Rotation Simulation Start" + " " * 17 + "#")
    print("#" + " " * 68 + "#")
    print("#" * 70 + "\n")

    input_structure = params["input_structure"]
    print(f"Reading input structure: {input_structure}")
    atoms = read(input_structure)

    print("\n" + "=" * 70)
    print("Step 1: Optimize the initial structure")
    print("=" * 70)
    atoms_opt, E0 = optimize_structure(atoms, params, silent=False)
    print(f"Initial optimized energy E0 = {E0:.6f} eV")
    print("Initial optimization files: initial_optimization.traj, initial_optimization.log")

    write("model.xyz", atoms_opt)

    with open("mc_energy.log", "w", encoding="utf-8") as f:
        f.write("# Step  Energy(eV)  DeltaE(eV)  Status\n")
        f.write(f"0  {E0:.6f}  0.000000  Initial\n")

    mc_stats = {
        "total_attempts": 0,
        "accepted_lower": 0,
        "accepted_metropolis": 0,
        "rejected": 0,
    }

    rotate_max = params["rotate_max"]
    temperature = params["temperature"]
    early_stop = params["early_stop"]

    print("\n" + "=" * 70)
    print("Monte Carlo Parameters")
    print("=" * 70)
    print(f"Max rotation steps:    {rotate_max}")
    print(f"Temperature:           {temperature} K")
    print(f"Early stop:            {early_stop}")
    print(f"Current reference E0:  {E0:.6f} eV")
    print("=" * 70)
    print("\nNote: optimization inside the MC loop does not write traj/log files.")
    print("=" * 70)

    for step in range(1, rotate_max + 1):
        print("\n" + "=" * 70)
        print(f"MC Step {step}/{rotate_max}")
        print("=" * 70)

        active_indices = find_active_carbon_atoms(atoms_opt, params, verbose=True)

        if len(active_indices) == 0:
            print("\nError: no active carbon atoms were found.")
            with open("mc_summary.log", "w", encoding="utf-8") as f:
                f.write(f"Monte Carlo simulation stopped at step {step}\n")
                f.write("Reason: No active carbon atoms found\n")
                f.write("\nStatistics:\n")
                f.write(f"Total attempts: {mc_stats['total_attempts']}\n")
                f.write(f"Accepted (lower energy): {mc_stats['accepted_lower']}\n")
                f.write(f"Accepted (Metropolis): {mc_stats['accepted_metropolis']}\n")
                f.write(f"Rejected: {mc_stats['rejected']}\n")
            raise RuntimeError("No active carbon atoms found. MC simulation terminated.")

        print("Applying trial rotation...")
        atoms_rotated, success = rotate_carbon_bond(atoms_opt, active_indices, params, verbose=True)

        if not success:
            print("Rotation failed. Skipping this MC step.")
            mc_stats["rejected"] += 1
            continue

        print("Optimizing the rotated structure in silent mode...")
        atoms_new, E_new = optimize_structure(atoms_rotated, params, silent=True)

        print(f"New energy E_new = {E_new:.6f} eV")
        delta_E = E_new - E0
        print(f"Energy difference dE = {delta_E:.6f} eV")

        mc_stats["total_attempts"] += 1

        accepted = False
        status = ""

        if E_new < E0:
            print("Accepted: lower energy.")
            accepted = True
            status = "Accepted(Lower)"
            mc_stats["accepted_lower"] += 1
        else:
            P_accept = np.exp(-delta_E / (kB * temperature))
            rand_num = random.random()

            print(f"Acceptance probability: P = {P_accept:.6f}")
            print(f"Random number:         r = {rand_num:.6f}")

            if rand_num < P_accept:
                print("Accepted by Metropolis criterion.")
                accepted = True
                status = "Accepted(Metropolis)"
                mc_stats["accepted_metropolis"] += 1
            else:
                print("Rejected.")
                status = "Rejected"
                mc_stats["rejected"] += 1

        with open("mc_energy.log", "a", encoding="utf-8") as f:
            f.write(f"{step}  {E_new:.6f}  {delta_E:.6f}  {status}\n")

        if accepted:
            atoms_opt = atoms_new
            E0 = E_new
            write("model.xyz", atoms_opt)
            print("Updated accepted structure written to model.xyz")
            print(f"Updated reference energy E0 = {E0:.6f} eV")

            if early_stop:
                print("\nEarly stop is enabled. MC loop terminates after the first accepted move.")
                break

        print("=" * 70)

    final_step = step
    print("\n" + "=" * 70)
    print("Monte Carlo Simulation Complete")
    print("=" * 70)

    with open("mc_summary.log", "w", encoding="utf-8") as f:
        f.write("Monte Carlo Rotation Simulation Summary\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total MC steps: {final_step}\n")
        f.write(f"Final energy: {E0:.6f} eV\n")
        f.write(f"Temperature: {temperature} K\n")
        f.write(f"Early stop: {early_stop}\n")
        f.write("\nStatistics:\n")
        f.write(f"Total attempts: {mc_stats['total_attempts']}\n")
        f.write(f"Accepted (lower energy): {mc_stats['accepted_lower']}\n")
        f.write(f"Accepted (Metropolis): {mc_stats['accepted_metropolis']}\n")
        f.write(f"Rejected: {mc_stats['rejected']}\n")
        if mc_stats["total_attempts"] > 0:
            accept_rate = (
                mc_stats["accepted_lower"] + mc_stats["accepted_metropolis"]
            ) / mc_stats["total_attempts"]
            f.write(f"Overall acceptance rate: {accept_rate:.2%}\n")

    print(f"Total MC steps:        {final_step}")
    print(f"Final energy:          {E0:.6f} eV")
    print(f"Total attempts:        {mc_stats['total_attempts']}")
    print(f"Accepted (lower):      {mc_stats['accepted_lower']}")
    print(f"Accepted (Metropolis): {mc_stats['accepted_metropolis']}")
    print(f"Rejected:              {mc_stats['rejected']}")

    if mc_stats["total_attempts"] > 0:
        accept_rate = (
            mc_stats["accepted_lower"] + mc_stats["accepted_metropolis"]
        ) / mc_stats["total_attempts"]
        print(f"Overall acceptance:    {accept_rate:.2%}")

    print("\nOutput files:")
    print("  - model.xyz                 : Latest accepted structure")
    print("  - mc_energy.log             : Monte Carlo energy history")
    print("  - mc_summary.log            : Simulation summary")
    print("  - initial_optimization.traj : Initial optimization trajectory")
    print("  - initial_optimization.log  : Initial optimization log")
    print("=" * 70)

    print("\n" + "#" * 70)
    print("#" + " " * 68 + "#")
    print("#" + " " * 19 + "Monte Carlo Rotation Simulation End" + " " * 18 + "#")
    print("#" + " " * 68 + "#")
    print("#" * 70 + "\n")

    return {
        "final_atoms": atoms_opt,
        "final_energy": E0,
        "mc_stats": mc_stats,
        "total_steps": final_step,
    }


def main():
    """
    Run the complete Monte Carlo rotation workflow.
    """
    params = read_edkmc_input("edkmc.in")
    results = monte_carlo_rotation(params)
    return results


if __name__ == "__main__":
    results = main()
