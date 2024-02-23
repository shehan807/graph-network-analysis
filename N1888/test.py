#!/usr/bin/env python
import networkx as nx
import mdtraj as md
import sys
import numpy as np
from joblib import Parallel, delayed
import multiprocessing as mp
import matplotlib.pyplot as plt
import pickle

def get_atoms(traj, res_list=["OCT"]):
    # TODO: Generalize to multi-residue
    atom_dict = {}
    for entry in res_list: 
        atom_dict[entry] = [[],[]]
        for res in traj.top.residues:
            if res.name == entry:
                atom_dict[entry][0].append(res.index)
                for atom in res.atoms:
                    atom_dict[entry][1].append(atom.index)
    return atom_dict

def get_xyz(traj, atoms):
    return traj.xyz[:, atoms, :]

def check_coords(xyz, box):
    # atom indices for positions larger than box length
    inds = np.where(xyz >= box)
    # move atoms back into box
    xyz[inds] -= box 
    # atom indices for positions less than box length 
    inds = np.where(xyz < 0.0)
    # move atoms back into box 
    xyz[inds] += box
    return xyz

def set_cell(ix, iy, iz, cell_size):
    """
    Parameters 
    -----------
    ix : int
        Id of cell along x direction
    iy : int
        Id of cell along y direction
    iz : int
        Id of cell along z direction
    cell_size : int
        Number of cells along each dimension

    Returns
    -----------
    cell_id : int
        cell id 
    """
    cell_id = 1 + (ix - 1 + cell_size) % cell_size + (iy - 1 + cell_size) % cell_size * cell_size + (iz - 1 + cell_size) % cell_size * cell_size**2
    return cell_id

def make_map(cell_size):
    """
    Parameters
    -----------
    cell_size : int
        Number of cells that each dimension is divided into

    Returns 
    -----------
    map : np.ndarray
        Contains int identifier for each cell
    """

    #Total number of cells 
    cell_num = cell_size**3

    #Neighboring cells for each cell
    map_size = 13 * cell_num

    #Initialize map
    cell_map = np.zeros((map_size + 1))

    #Loop through the number of cells in the x, y and z dimensions. Do this with Fortran indexing (starting from 1) to match old code
    for iz in range(1, cell_size+1):
        for iy in range(1, cell_size+1):
            for ix in range(1, cell_size+1):
                #Cell id
                cell_id = (set_cell(ix, iy, iz, cell_size) - 1) * 13
                #Set the next 13 entries based off of the cell id
                cell_map[cell_id+1] = set_cell(ix+1, iy, iz, cell_size)
                cell_map[cell_id+2] = set_cell(ix + 1, iy + 1, iz, cell_size)
                cell_map[cell_id+3] = set_cell(ix, iy+1, iz, cell_size)
                cell_map[cell_id+4] = set_cell(ix - 1, iy+1, iz, cell_size)
                cell_map[cell_id+5] = set_cell(ix + 1, iy, iz - 1 , cell_size)
                cell_map[cell_id+6] = set_cell(ix + 1, iy + 1, iz - 1 , cell_size)
                cell_map[cell_id+7] = set_cell(ix, iy + 1, iz - 1, cell_size)
                cell_map[cell_id+8] = set_cell(ix - 1, iy + 1, iz - 1 , cell_size)
                cell_map[cell_id+9] = set_cell(ix + 1, iy, iz + 1 , cell_size)
                cell_map[cell_id+10] = set_cell(ix + 1, iy + 1, iz + 1 , cell_size)
                cell_map[cell_id+11] = set_cell(ix, iy + 1, iz + 1 , cell_size)
                cell_map[cell_id+12] = set_cell(ix - 1, iy + 1, iz + 1 , cell_size)
                cell_map[cell_id+13] = set_cell(ix, iy, iz + 1 , cell_size)

    return cell_map

def make_head(head_list, o_list, xyz, box, cell_size, cutoff):

    # dimension of each cell 
    cl = box / cell_size

    # cell size must be larger than cutoff (i.e., ~coordination distance)
    if cl < cutoff:
        print("Error, cell size too small for cutoff!!")

    for i_atom, i in enumerate(range(0, xyz.shape[0], 3)):
        icell = 1 + int(xyz[i,0] / cl)\
                  + int(xyz[i,1] / cl) * cell_size\
                  + int(xyz[i,2] / cl) * cell_size**2\

        o_list[i_atom+1] = head_list[icell]
        head_list[icell] = i_atom + 1
        
        return head_list, o_list


def get_network(xyz, box, cell_map, cell_size, cutoff, frame):
    print(frame)

    # total number of cells
    cell_num = cell_size**3

    # head of cell linked list
    head_list = np.zeros(cell_num+1)

    # link each atom through the head list in each cell 
    o_list = np.zeros(int(xyz.shape[0] / 3) + 1)

    # fill head and o_list
    head_list, o_list = make_head(head_list, o_list, xyz, box, cell_size, cutoff)

    # determine h bonds/edges 
    edges = get_frame_edges(xyz, box, cell_map, cell_num, head_list, o_list)

    return edges

# Processes count
cores = 12 
# Load trj to MDTraj
traj = md.load("npt_final-28.dcd", top="npt_final-28.pdb") 
# residue and atom index list for relevant groups
atom_dict = get_atoms(traj, ["OCT"])
cation_residues, cation_atoms = atom_dict["OCT"][0], atom_dict["OCT"][0]
# index list from 0
res_index = [i for i in range(len(cation_residues))]

frames = len(traj)

# np array of cation atoms
cation_xyz = get_xyz(traj, cation_atoms)
# check for out-of-box coords and remove
cation_xyz = check_coords(cation_xyz, traj.unitcell_lengths[0,0])
# divide box into cells 
cell_map = make_map(10)

num_frames = len(traj)





