import mdtraj as md
import os, sys
import numpy as np
import warnings
import pickle 

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
        traj_filtered = traj.atom_slice(atoms)
        with open(trj_fil_file, "wb") as file: 
            pickle.dump(traj_filtered, file)

    res_atoms = list(traj_filtered.topology.atoms)
    atom_per_res = int(len(atoms) / len(residues))
    xyz = get_xyz(traj, atoms, box, check_coords=True)
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

def get_angle(xyz, i_atom, j_atom, h_atom, box):
    """
    Parameters
    -----------
    xyz : np.ndarray
        contains the coordinates of water molecules
    i_oxygen : int
        index of oxygen atom on molecule i
    j_oxygen : int
        index of oxygen atom on molecule j
    h_atom : int
        index of hydrogen atom that the angle is computed for
    box : float
        box dimension

    Returns
    -----------
    theta : float
        computed angle
    """
    #Displacement between i_oxygen and the h_atom
    dist, disp_1 = get_distance(water_xyz[i_oxygen, :], water_xyz[h_atom, :], box)
    #Displacement between i_oxygen and j_oxygen
    dist, disp_2 = get_distance(water_xyz[i_oxygen, :], water_xyz[j_oxygen, :], box)

    #Compute cosine as dot product
    mag_u = np.sqrt((disp_1 * disp_1).sum())
    mag_v = np.sqrt((disp_2 * disp_2).sum())
    cosine = ((disp_1 * disp_2).sum()) / mag_u / mag_v
    
    if cosine < -0.999999999:
      theta = 3.14159265359
    elif cosine > 0.999999999:
      theta = 0.0
    else:
      theta = np.arccos(cosine)
    theta = theta*180/3.14159265359

    return theta

def check_hbond_criteria(water_xyz, i_oxygen, j_oxygen, box):
    """
    Parameters
    -----------
    water_xyz : np.ndarray
        contains the coordinates of water molecules
    i_oxygen : int
        index of oxygen atom on molecule i
    j_oxygen : int
        index of oxygen atom on molecule j
    box : float
        box dimension

    Returns
    -----------
    bool
        Returns True if the hydrogen bond criteria is met, False otherwise
    """

    #Hydrogens on i_oxygen
    h1 = i_oxygen + 1
    h2 = i_oxygen + 2
    h_atoms = np.asarray([h1, h2])

    #Get distances between hydrogens and j_oxygen
    dist_1, disp = get_distance(water_xyz[h1, :], water_xyz[j_oxygen, :], box)
    dist_2, disp = get_distance(water_xyz[h2, :], water_xyz[j_oxygen, :], box)
    dists = np.asarray([dist_1, dist_2])

    #Find the minimum H-O distance
    h_min = np.argmin(dists)
    dist = dists[h_min]

    #Check criteria
    if dist < 0.245:
        #Check angle criteria
        angle = get_angle(water_xyz, i_oxygen, j_oxygen, h_atoms[h_min], box)
        if angle < 30:
            return True

    #Hydrogens on j_oxygen
    h1 = j_oxygen + 1
    h2 = j_oxygen + 2
    h_atoms = np.asarray([h1, h2])

    #Get distances between hydrogens and i_oxygen
    dist_1, disp = get_distance(water_xyz[h1, :], water_xyz[i_oxygen, :], box)
    dist_2, disp = get_distance(water_xyz[h2, :], water_xyz[i_oxygen, :], box)
    dists = np.asarray([dist_1, dist_2])

    #Find the minimum H-O distance
    h_min = np.argmin(dists)
    dist = dists[h_min]
    
    #Check criteria
    if dist < 0.245:
        #Check angle criteria
        angle = get_angle(water_xyz, j_oxygen, i_oxygen, h_atoms[h_min], box)
        if angle < 30:
            return True

    return False

def check_criteria(criteria, filtered_atoms, i_mol, i_mol_atom_ind, i_mol_atoms, j_mol, j_mol_atom_ind, j_mol_atoms, xyz, traj_filtered):
    """
    Parameters
    ----------
    criteria : dict
    filtered_atoms : list
    i_mol : int
    i_mol_atoms : list
    j_mol : int
    j_mol_atoms : list
    xyz : list

    Returns: 
    --------
    chk : bool
        boolean for whether criteria is satisfied or not
    """

    chk = False # NOTE: this could be a source of some issue
    dist_pairs = []
    for criterion in criteria:
        if criterion["distance"] is not None: 
            # get names of residue IDs and atoms being evaluated 
            i_atom_name, j_atom_name = criterion["name"].split('-')
            i_atom = int(i_mol_atom_ind + i_mol_atoms.index(i_atom_name))
            j_atom = int(j_mol_atom_ind + j_mol_atoms.index(j_atom_name))
            
            dist_pairs.append([i_atom, j_atom])
            # obtain box dimension
            # box = traj_filtered.unitcell_lengths[0,0]
    #print("dist_pairs")
    #print(dist_pairs)

    # get distance between i_atom, j_atom
    dist = md.compute_distances(traj_filtered, dist_pairs) 
    #print(dist)
    
    # check criteria
    dist_cri_chks = np.sum(dist < criterion["distance"])
    chk_dist = dist_cri_chks >= criterion["min_true"] 
    #print("checking criterion")
    #print(dist < criterion["distance"])
    #print(chk_dist)
    # once an atom pair satisfies criteria, the molecules form an edge
    if chk_dist:
        #print(f"Criteria met!")
        #print(f"Criteria met {traj_filtered.topology.atom(i_atom)}---{traj_filtered.topology.atom(j_atom)} dist = {dist}!!!")
        #print(dist_pairs)
        #print(dist)
        chk = chk_dist or chk
    return chk
