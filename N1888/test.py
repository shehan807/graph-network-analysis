#!/usr/bin/env python
import networkx as nx
import mdtraj as md
import sys
import numpy as np
from joblib import Parallel, delayed
import multiprocessing as mp
import matplotlib.pyplot as plt
import pickle

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
    #print(len(atom_dict[entry][0]))
    #print(len(atom_dict[entry][1]))
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
                cell_map[cell_id +  1] = set_cell(ix + 1, iy    , iz    , cell_size)
                cell_map[cell_id +  2] = set_cell(ix + 1, iy + 1, iz    , cell_size)
                cell_map[cell_id +  3] = set_cell(ix    , iy + 1, iz    , cell_size)
                cell_map[cell_id +  4] = set_cell(ix - 1, iy + 1, iz    , cell_size)
                cell_map[cell_id +  5] = set_cell(ix + 1, iy    , iz - 1, cell_size)
                cell_map[cell_id +  6] = set_cell(ix + 1, iy + 1, iz - 1, cell_size)
                cell_map[cell_id +  7] = set_cell(ix    , iy + 1, iz - 1, cell_size)
                cell_map[cell_id +  8] = set_cell(ix - 1, iy + 1, iz - 1, cell_size)
                cell_map[cell_id +  9] = set_cell(ix + 1, iy    , iz + 1, cell_size)
                cell_map[cell_id + 10] = set_cell(ix + 1, iy + 1, iz + 1, cell_size)
                cell_map[cell_id + 11] = set_cell(ix    , iy + 1, iz + 1, cell_size)
                cell_map[cell_id + 12] = set_cell(ix - 1, iy + 1, iz + 1, cell_size)
                cell_map[cell_id + 13] = set_cell(ix    , iy    , iz + 1, cell_size)

    return cell_map

def make_head(head_list, o_list, xyz, box, cell_size, cutoff):

    # dimension of each cell 
    cl = box / cell_size

    # cell size must be larger than cutoff (i.e., ~coordination distance)
    if cl < cutoff:
        print("Error, cell size too small for cutoff!!")
    for i_atom, i in enumerate(range(0, xyz.shape[0], 106)):
        icell = 1 + int(xyz[i,0] / cl)\
                  + int(xyz[i,1] / cl) * cell_size\
                  + int(xyz[i,2] / cl) * cell_size**2
        o_list[i_atom] = head_list[icell]
        head_list[icell] = i_atom + 1
    return head_list, o_list

def get_frame_edges(xyz, box, cell_map, cell_num, head_list, o_list):
    """
    Parameters
    -----------
    xyz : np.ndarray
        contains the positions of all relevant atoms
    box : float
        box dimension
    cell_map : np.ndarray
        Contains the divided map of the cell
    cell_num : int
        Total number of cells
    head_list : np.ndarray
        The head list of molecules in each cell
    o_list : np.ndarray
        Linked list for each molecule

    Returns
    -----------
    edges : list
        Contains pairs of 'connected' molecules
    """
    edges = []
    for i in range(1, cell_num+1):
        #Get initial molecule from the head list
        i_mol = head_list[i]
        #Some cells won't have a molecule, so only need to loop through the cell if i_mol isn't 0
        while (i_mol > 0):
            #Get the linked j_atom from o_list
            j_mol = o_list[int(i_mol)]
            #print(i_mol, j_mol)

            #If j_mol is nonzero
            while (j_mol > 0):
                for C1 in range(0,26):
                    for C2 in range(0,26):
                        i_carbon = int((i_mol - 1)*C1)
                        j_carbon = int((j_mol - 1)*C2)
                        #Get distance of this pair
                        dist, disp = get_distance(xyz[i_carbon], xyz[j_carbon], box)
                        #Check O-O criteria
                        if dist < 0.70 and dist > 0.30:
                            if (int(i_mol)-1, int(j_mol)-1) not in edges:
                                edges.append((int(i_mol)-1, int(j_mol)-1))
                                break
                        #print(f"C{C}-C{C} Dist = {dist} nm")

                #See if there's a molecule linked to j_mol
                j_mol = o_list[int(j_mol)]

            #Now loop through neighboring cells and see if there are any molecules bonded to i_mol
            jcell0 = 13 * (i - 1)
            for neighbor in range(1, 14):
                #Get cell number for neighboring cell of cell i
                jcell = int(cell_map[jcell0 + neighbor])
                #Head molecule index of jcell
                j_mol = head_list[jcell]

                while (j_mol > 0):
                    
                    for C1 in range(0,26):
                        for C2 in range(0,26):
                            i_carbon = int((i_mol - 1)*C1)
                            j_carbon = int((j_mol - 1)*C2)
                            dist, disp = get_distance(xyz[i_carbon], xyz[j_carbon], box)
                            #Check O-O criteria (must be less than 0.36 nm)
                            if dist < 0.70 and dist > 0.30:
                                if (int(i_mol)-1, int(j_mol)-1) not in edges:
                                    edges.append((int(i_mol)-1, int(j_mol)-1))
                                    break
                            #print(f"C{C}-C{C} Dist = {dist} nm")
                            #break
                    #See if there's an atom linked to j_mol
                    j_mol = o_list[int(j_mol)]
            
            #See if there's a linked molecule to i_mol
            i_mol = o_list[int(i_mol)]

    return edges


def get_network(xyz, box, cell_map, cell_size, cutoff, frame):
    print(frame)

    # total number of cells
    cell_num = cell_size**3

    # head of cell linked list
    head_list = np.zeros(cell_num+1)
    #print(cell_num)
    #print(len(head_list))
    # link each atom through the head list in each cell 
    o_list = np.zeros(int(xyz.shape[0] / 3) + 1)

    # fill head and o_list
    head_list, o_list = make_head(head_list, o_list, xyz, box, cell_size, cutoff)

    # determine h bonds/edges 
    edges = get_frame_edges(xyz, box, cell_map, cell_num, head_list, o_list)

    return edges

def make_graph(edge, res_index):
    """
    Parameters
    -----------
    edge : list
        contains hydrogen bond pairs
    res_index : list
        index of each water molecule

    Returns
    -----------
    graph : NetworkX Graph
        graph object
    """
    graph = nx.Graph()
    for i in res_index: graph.add_node(i)
    graph.add_edges_from(edge)
    return graph

def compute_metric(graph):
    """
    Parameters
    -----------
    graph : Networkx graph

    Returns
    -----------
    diam_water : list
        List of each water cluster diameter
    """

    #Compute diameter of the cluster here
    diam_water = []
    for cl in nx.connected_components(graph):
        diam_water.append(nx.diameter(graph.subgraph(cl)))
    return diam_water

def plt_metric(metrics):
    
    diams = []
    for metric in metrics:
        for m in metric: diams.append(m)

    weights_water = np.ones_like(diams)/float(len(diams))
    water_hist, water_edges, water_p = plt.hist(diams, weights=weights_water, bins=40, alpha=0.5, color='midnightblue',label='N1888+')
    plt.legend()
    plt.xlabel("Diameter Length")
    plt.ylabel("Probability")

    plt.savefig("diam_n1888-combinations.png")

