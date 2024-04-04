import mdtraj as md
import sys
import numpy as np
import warnings

def initialize(dcd, pdb, residue_name, num_cells, cutoff):
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
    atom_per_res = int(len(atoms) / len(residues))
    xyz = get_xyz(traj, atoms, box, check_coords=True)
    print(len(xyz))
    #print(len(xyz[0])/len(residues))
    
    return xyz, traj, residues, atoms, num_frames, box, atom_per_res
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
    #print(f"({len(atoms)/len(residues)}) atoms in residue")
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
        print(f"Found {len(inds_g)} positions greater than box dim.")
        print(f"Found {len(inds_l)} positions less than box dim.")
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
    #print(f"xyz.shape = {xyz.shape}.")
    linked_list = np.zeros(int(num_atoms / atom_per_res))
    print(atom_per_res)
    #print(f"cell length = {cell_length} (cutoff = {cutoff}).")
    #print(f"head list length = {len(head_list)}")
    #print(f"linked list length = {len(linked_list)}")
    if cell_length < cutoff: # cell size must be larger than cutoff (i.e., coordination distance)
        warnings.warn("Cell size too small for cutoff!!", UserWarning)
    
    for iatom, i in enumerate(range(0, num_atoms, atom_per_res)):
        icell = int(xyz[i, 0] / cell_length)\
              + int(xyz[i, 1] / cell_length) * num_cells\
              + int(xyz[i, 2] / cell_length) * num_cells**2
        #print(f"atom index {iatom} located in cell {icell}")
        linked_list[iatom] = head_list[icell]
        #print(f"linked_list[{iatom}] = head_list[{icell}] = {head_list[icell]}.")
        head_list[icell] = iatom
        #print(f"head_list[{icell}] = {iatom}")


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
    print(f"Creating cell_map of {len(cell_map)} elements, i.e., 13 for each of the {total_cells} cells.")
    offsets = [(1, 0, 0), (1, 1, 0), (0, 1, 0), (-1, 1, 0),
               (1, 0, -1), (1, 1, -1), (0, 1, -1), (-1, 1, -1),
               (1, 0, 1), (1, 1, 1), (0, 1, 1), (-1, 1, 1),
               (0, 0, 1)
               ]
    for iz in range(0,num_cells+0):
        for iy in range(0,num_cells+0):
            for ix in range(0,num_cells+0):
                icell = set_cell(ix, iy, iz, num_cells) * 13
                # print(f"%%% ({ix},{iy},{iz}) icell = {icell} %%%")
                for map_index, (dx, dy, dz) in enumerate(offsets): 
                    neighbor_index = set_cell(ix + dx, iy + dy, iz + dz, num_cells)
                    cell_map[icell + map_index + 0] = neighbor_index
                    #print(f"({ix+dx},{iy+dy},{iz+dz}) icell = {neighbor_index}")
                    #print(f"cell_map[{icell+map_index+0}]={neighbor_index}")
    return cell_map 


def get_distance():
    pass
def get_angle():
    pass
def hbond_criteria():
    pass

