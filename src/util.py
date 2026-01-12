import mdtraj as md
import os, sys
import numpy as np
import warnings
import pickle 
from copy import deepcopy 
import logging

def initialize(dcd, pdb, out, residue_name, num_cells, cutoff):
    """Initializes parameters for networkx functions.
    Parameters
    ----------
    dcd : 
        .dcd trajectory file
    pdb : 
        .pdb topology file
    residue_names : 
        species name to match in pdb file for inter-species analysis
    
    Returns
    -------
    xyz : 
        ...
    traj : 
        ...
    cell_map : 
        ...
    residues : 
        ...
    atoms : 
        ...
    num_frames : 
        ...
    """
    traj = md.load(dcd, top=pdb)
    num_frames = len(traj) 
    box = traj.unitcell_lengths[0,0]
    residues, atoms = get_atoms(traj, residue_name)
   
    # limit traj object slicing 
    trj_fil_file = os.path.join(out, "trj_filtered.pkl")
    if os.path.isfile(trj_fil_file):
        with open(trj_fil_file, "rb") as file:
            traj_filtered = pickle.load(file)
    else: 
        logging.debug("No filtered trajectory pkl file exists, creating one now.")
        traj_filtered = traj.atom_slice(atoms)
        with open(trj_fil_file, "wb") as file: 
            pickle.dump(traj_filtered, file)

    res_atoms = list(traj_filtered.topology.atoms)
    atom_per_res = int(len(atoms) / len(residues))
    xyz = get_xyz(traj, atoms, box, check_coords=True)
    logging.debug("Produced trajectory length of %d for %d atoms.",len(xyz), len(atoms))

    return xyz, traj, traj_filtered, residues, atoms, num_frames, box, atom_per_res

def get_atoms(traj, residue_name):
    """Filter out list of indices for residue and respective atoms
    Parameters
    ----------
    traj : mdtraj.core.trajectory.Trajectory 
        MDTraj trajectory object
    residue_name : str
        name of residue of interest matching pdb topology file 
    Returns
    -------
    residues : list  
        list of residue indices for selected residue name
    atoms : list 
        list of atom indices for selected atom name 
    
    """
    residues, atoms = [], []
    for res in traj.top.residues: 
        if res.name == residue_name: 
            residues.append(res.index)
            for atom in res.atoms: 
                atoms.append(atom.index)
    return residues, atoms

def get_xyz(traj, atoms, box, check_coords=True):
    """
    Parameters
    ----------
    traj : 
    atoms : 
    box : 
    check_coords : 
    Returns
    -------
    xyz : 
    """
    xyz = traj.xyz[:, atoms, :]
    if check_coords: 
        inds_g = np.where(xyz >= box)
        inds_l = np.where(xyz < 0.0)
        xyz[inds_g] -= box 
        xyz[inds_l] += box
    return xyz

def make_head(xyz, box, num_cells, cutoff, atom_per_res):
    """
    Returns
    -------
    head_list : 
        `head-of-chain` array, which has one element for each cell containing the identification number of one of the molecules sorted into that cell
    linked_list : 
        `linked-list` array1, an list where the identification number (from the head list) is used to address an element of linked_list which contains the number of the next molecules in that cell, i.e., linked_list elements for a specified molecule is the indec of the next molecule in that cell, etc.

    """
    total_cells = num_cells*num_cells*num_cells
    cell_length = box / num_cells
    num_atoms = xyz.shape[0]

    head_list = np.zeros(total_cells)
    linked_list = np.zeros(int(num_atoms / atom_per_res))
    if cell_length < cutoff: # cell size must be larger than cutoff (i.e., coordination distance)
        warnings.warn("Cell size too small for cutoff!!", UserWarning)
    
    for iatom, i in enumerate(range(0, num_atoms, atom_per_res)):
        icell = int(xyz[i, 0] / cell_length)\
              + int(xyz[i, 1] / cell_length) * num_cells\
              + int(xyz[i, 2] / cell_length) * num_cells**2
        linked_list[iatom] = head_list[icell]
        head_list[icell] = iatom
    return head_list, linked_list

def set_cell(ix, iy, iz, num_cells):
    icell = (ix + num_cells) % num_cells\
          + (iy + num_cells) % num_cells * num_cells\
          + (iz + num_cells) % num_cells * num_cells**2
    return icell

def make_map(num_cells):
    """...
    Creates cell_map, a list of 13 neighboring cells of each of the small cells in the central box including periodic boundary conditions (i.e., to get all 26 neighbors of a cube).
    See Reference: https://github.com/glennklockwood/allen-tildesley/blob/master/link-cell.f
    Parameters
    ----------
    num_cells : int
        number of cells each box dimension is divided into to create an MxMxM lattice of cells
    Returns
    -------
    map : np.ndarray
        ...
    """
    total_cells = num_cells*num_cells*num_cells
    map_size = total_cells * 13 
    cell_map = np.zeros(map_size + 0) # exclude (0,0,0) 
    offsets = [(1, 0, 0), (1, 1, 0), (0, 1, 0), (-1, 1, 0),
               (1, 0, -1), (1, 1, -1), (0, 1, -1), (-1, 1, -1),
               (1, 0, 1), (1, 1, 1), (0, 1, 1), (-1, 1, 1),
               (0, 0, 1)
               ]
    for iz in range(0,num_cells+0):
        for iy in range(0,num_cells+0):
            for ix in range(0,num_cells+0):
                icell = set_cell(ix, iy, iz, num_cells) * 13
                for map_index, (dx, dy, dz) in enumerate(offsets): 
                    neighbor_index = set_cell(ix + dx, iy + dy, iz + dz, num_cells)
                    cell_map[icell + map_index + 0] = neighbor_index
    return cell_map 

def get_distance(i_atom, j_atom, box):
    """
    Parameters 
    -----------
    i_atom : int
        index of atom i
    j_atom : int
        index of atom j
    box : float
        box dimension

    Returns
    -----------
    dist : float
        minimum image distance
    disp : float
        displacement vector (minimum imaged as well)
    """
    disp = i_atom - j_atom
    #See if the vector is more than half the box length
    shift = box * np.round(disp/box)
    disp -= shift
    return np.sqrt((disp*disp).sum()), disp

def compute_angle(bond_atom_xyz, vertex_atom_xyz, ref_atom_xyz, box):
    """
    Compute angle at vertex atom with periodic boundary conditions.

    Parameters
    ----------
    bond_atom_xyz : np.ndarray
        XYZ coordinates of the bond atom on i_mol (shape: (3,))
    vertex_atom_xyz : np.ndarray
        XYZ coordinates of the vertex atom where angle is measured (shape: (3,))
    ref_atom_xyz : np.ndarray
        XYZ coordinates of the reference atom on j_mol (shape: (3,))
    box : float
        Cubic box dimension for periodic boundary conditions (nm)

    Returns
    -------
    angle : float
        Angle at vertex in degrees (range: 0-180)
        180 degrees = perfectly linear arrangement

    Notes
    -----
    Computes the angle: bond_atom --- vertex_atom --- ref_atom

    The angle is measured at vertex_atom, between:
      - Vector from vertex_atom to bond_atom (the intramolecular bond arm)
      - Vector from vertex_atom to ref_atom (the intermolecular reference arm)

    Uses minimum image convention for periodic boundaries.

    Example:
        For hydrogen bonding O-H...O where H is vertex:
        - bond_atom = O on donor molecule (bonded to H)
        - vertex_atom = H on donor molecule (the vertex)
        - ref_atom = O on acceptor molecule
        - An angle of 180 means O-H-O are collinear (ideal H-bond geometry)
    """
    # Compute displacement vector from vertex to bond atom (with PBC)
    _, vec_to_bond = get_distance(bond_atom_xyz, vertex_atom_xyz, box)

    # Compute displacement vector from vertex to reference atom (with PBC)
    _, vec_to_ref = get_distance(ref_atom_xyz, vertex_atom_xyz, box)

    # Compute magnitudes
    mag_bond = np.sqrt(np.sum(vec_to_bond * vec_to_bond))
    mag_ref = np.sqrt(np.sum(vec_to_ref * vec_to_ref))

    # Compute cosine of angle via dot product
    dot_product = np.sum(vec_to_bond * vec_to_ref)
    cosine = dot_product / (mag_bond * mag_ref)

    # Clamp cosine to [-1, 1] to handle numerical precision issues
    cosine = np.clip(cosine, -1.0, 1.0)

    # Convert to degrees
    angle = np.degrees(np.arccos(cosine))

    return angle

def check_criteria(criteria, filtered_atoms, i_mol, i_mol_atom_ind, i_mol_atoms,
                   j_mol, j_mol_atom_ind, j_mol_atoms, xyz, traj_filtered, box):
    """
    Check if molecule pair (i_mol, j_mol) satisfies edge criteria (with angle support).

    Parameters
    ----------
    criteria : list of dict
        List of criteria from config.yaml. Each criterion has:
        - name: "atom1-atom2" (atom1 on i_mol, atom2 on j_mol)
        - distance: float threshold in nm
        - min_true: int minimum passing pairs
        - angle: optional dict with {bond_atom, min, max} or None
    filtered_atoms : list
        List of atom names from filtered trajectory
    i_mol : int
        Index of molecule i
    i_mol_atom_ind : int
        Starting atom index for molecule i in filtered trajectory
    i_mol_atoms : list
        List of atom names in molecule i
    j_mol : int
        Index of molecule j
    j_mol_atom_ind : int
        Starting atom index for molecule j in filtered trajectory
    j_mol_atoms : list
        List of atom names in molecule j
    xyz : np.ndarray
        Coordinate array for current frame
    traj_filtered : mdtraj.Trajectory
        Filtered trajectory object for distance computation
    box : float
        Box dimension for periodic boundary conditions

    Returns
    -------
    chk : bool
        True if criteria are satisfied (edge should be created), False otherwise

    Notes
    -----
    For each criterion:
    - If only distance is specified: pair passes if distance < threshold
    - If distance + angle is specified: pair passes if BOTH distance < threshold
      AND angle is within [min, max] range

    Angle definition for criterion name "H1-O" with bond_atom "O":
    - Vertex: H1 on i_mol (first atom in name)
    - Bond arm: O on i_mol (specified by bond_atom)
    - Reference arm: O on j_mol (second atom in name)
    - Computed angle: O(i_mol) --- H1(i_mol) --- O(j_mol)
    - 180 degrees = perfectly linear (ideal H-bond geometry)

    Config example:
        criteria:
          - name: "H1-O"
            distance: 0.245
            min_true: 1
            angle:
              bond_atom: "O"
              min: 150.0
              max: 180.0
    """
    chk = False
    dist_pairs = []
    dist_pairs_names = []
    criteria_info = []  # store criterion info for each pair

    for criterion in criteria:
        if criterion["distance"] is not None:
            # get names of residue IDs and atoms being evaluated
            i_atom_name, j_atom_name = criterion["name"].split('-')
            i_atom = int(i_mol_atom_ind + i_mol_atoms.index(i_atom_name))
            j_atom = int(j_mol_atom_ind + j_mol_atoms.index(j_atom_name))

            dist_pairs.append([i_atom, j_atom])
            dist_pairs_names.append([traj_filtered.topology.atom(int(i_atom)),
                                     traj_filtered.topology.atom(int(j_atom))])
            criteria_info.append({
                'criterion': criterion,
                'i_atom': i_atom,
                'j_atom': j_atom,
                'i_atom_name': i_atom_name,
                'j_atom_name': j_atom_name
            })

    # return early if no distance pairs to check
    if len(dist_pairs) == 0:
        return chk

    # get distance between i_atom, j_atom
    dist = md.compute_distances(traj_filtered, dist_pairs)

    # check criteria with angle support
    passing_count = 0
    frame_xyz = traj_filtered.xyz[0]  # coordinates for angle computation

    for idx, info in enumerate(criteria_info):
        criterion = info['criterion']
        pair_dist = dist[0, idx]

        # check distance criterion
        if pair_dist >= criterion["distance"]:
            continue  # distance check failed

        # distance check passed - now check angle if configured
        angle_config = criterion.get("angle")

        # Handle angle: None, angle: "None", or missing angle field
        if angle_config is None or angle_config == "None" or not isinstance(angle_config, dict):
            # no angle criterion - distance-only check passed
            passing_count += 1
            logging.debug("Distance criterion passed: %s-%s dist=%.4f < %.4f",
                         info['i_atom_name'], info['j_atom_name'],
                         pair_dist, criterion["distance"])
        else:
            # angle criterion is configured
            bond_atom_name = angle_config.get("bond_atom")
            angle_min = angle_config.get("min", 0.0)
            angle_max = angle_config.get("max", 180.0)

            if bond_atom_name is None:
                logging.warning("Angle config missing 'bond_atom' for criterion %s",
                               criterion["name"])
                continue

            # find bond atom index on i_mol
            try:
                bond_atom_local = i_mol_atoms.index(bond_atom_name)
            except ValueError:
                logging.warning("Bond atom '%s' not found in i_mol atoms: %s",
                               bond_atom_name, i_mol_atoms)
                continue

            bond_atom = int(i_mol_atom_ind + bond_atom_local)

            # get coordinates for angle computation
            bond_atom_xyz = frame_xyz[bond_atom]
            vertex_atom_xyz = frame_xyz[info['i_atom']]  # i_atom is the vertex
            ref_atom_xyz = frame_xyz[info['j_atom']]

            # compute angle
            angle = compute_angle(bond_atom_xyz, vertex_atom_xyz, ref_atom_xyz, box)

            # check angle criterion
            if angle_min <= angle <= angle_max:
                passing_count += 1
                logging.debug("Distance+Angle criterion passed: %s-%s dist=%.4f angle=%.1f",
                             info['i_atom_name'], info['j_atom_name'], pair_dist, angle)
            else:
                logging.debug("Angle criterion failed: %s-%s angle=%.1f not in [%.1f, %.1f]",
                             info['i_atom_name'], info['j_atom_name'], angle, angle_min, angle_max)

    # check if enough pairs passed (using last criterion's min_true for backward compat)
    min_true = criteria_info[-1]['criterion'].get("min_true", 1) if criteria_info else 1
    chk_dist = passing_count >= min_true

    if passing_count > 0:
        logging.debug("%d atom pairs passed criteria between %s---%s",
                     passing_count, traj_filtered.topology.residue(int(i_mol)),
                     traj_filtered.topology.residue(int(j_mol)))

    # once enough atom pairs satisfy criteria, the molecules form an edge
    if chk_dist:
        chk = True
    return chk
