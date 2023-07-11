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

def get_angle(water_xyz, i_oxygen, j_oxygen, h_atom, box):
    """
    Parameters
    -----------
    water_xyz : np.ndarray
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

def get_frame_edges(water_xyz, box, cell_map, cell_num, head_list, o_list):
    """
    Parameters
    -----------
    water_xyz : np.ndarray
        contains the positions of all water atoms
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
        Contains pairs of hydrogen bonded water molecules
    """
    print(" %%% printing head_list %%%")
    print(head_list)
    print(" %%% printing cell_num %%%")
    print(cell_num)
    print(" %%% printing head_list[1] %%%")
    print(head_list[1])
    edges = []
    for i in range(1, cell_num+1):
        #Get initial molecule from the head list
        i_mol = head_list[i]

        #Some cells won't have a molecule, so only need to loop through the cell if i_mol isn't 0
        while (i_mol > 0):
            #Get the linked j_atom from o_list
            j_mol = o_list[int(i_mol)]

            #If j_mol is nonzero
            while (j_mol > 0):
                #Oxygen index of i_mol
                i_oxygen = int((i_mol - 1)*3)
                #Oxygen index of j_mol
                j_oxygen = int((j_mol - 1)*3)
                #Get distance of this pair
                dist, disp = get_distance(water_xyz[i_oxygen], water_xyz[j_oxygen], box)

                #Check O-O criteria
                if dist < 0.36:
                    #Check O-H and O-H-O criteria to determine if the molecules are hydrogen bonded
                    if check_hbond_criteria(water_xyz, i_oxygen, j_oxygen, box): edges.append((int(i_mol)-1, int(j_mol)-1))

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
                    i_oxygen = int((i_mol - 1)*3)
                    j_oxygen = int((j_mol - 1)*3)
                    dist, disp = get_distance(water_xyz[i_oxygen], water_xyz[j_oxygen], box)

                    #Check O-O criteria (must be less than 0.36 nm)
                    if dist < 0.36:
                        #Check H-O (< 0.245 nm) and O-H-O criteria (less than 30 degrees)
                        if check_hbond_criteria(water_xyz, i_oxygen, j_oxygen, box): edges.append((int(i_mol)-1, int(j_mol)-1))

                    #See if there's an atom linked to j_mol
                    j_mol = o_list[int(j_mol)]
            
            #See if there's a linked molecule to i_mol
            i_mol = o_list[int(i_mol)]

    return edges

def make_head(head_list, o_list, water_xyz, box, cell_size, cutoff):
    """
    Parameters
    -----------
    head_list : np.ndarray
        head of cell linked list
    o_list : np.ndarray
        contains the head value for each atom
    water_xyz : np.ndarray
        contains the positions of all water atoms
    box : float
        box dimension
    cell_size : int
        Number of cells along each dimension
    cutoff : float
        The box dimension divided by the cell_size must be larger than this amount. Don't adjust this as it should probably be larger than the O-O hydrogen bond distance

    Returns
    -----------
    head_list : np.ndarray
        Filled head_list
    o_list : np.ndarray
        Filled o_list
    """
    #Dimension in each direction for the cells in the box
    cl = box/cell_size

    #The dimension on the side of each cell should be larger than the cutoff value
    if cl < cutoff:
        print("Error, cell size too small for cutoff")

    #Loop through the number of water molecules, compute the cell that each molecule oxygen is in and fill in the o_list and head_list
    for i_atom, i in enumerate(range(0, water_xyz.shape[0], 3)):
        icell = 1 + int(water_xyz[i, 0]/cl) + int(water_xyz[i, 1]/cl) * cell_size + int(water_xyz[i, 2]/cl) * cell_size**2
        o_list[i_atom+1] = head_list[icell]
        head_list[icell] = i_atom + 1

    return head_list, o_list

def get_network(water_xyz, box, cell_map, cell_size, cutoff, frame):
    """
    Parameters
    -----------
    water_xyz : np.ndarray
        Array containing only water positions
    box : float
        Box dimension
    cell_map : np.ndarray
        Contains the divided map of the cell
    cell_size : int
        Number of cells along each dimension. Adjustable depending on how many cells you want
    cutoff : float
        The box dimension divided by the cell_size must be larger than this amount. Don't adjust this as it should probably be larger than the O-O hydrogen bond distance
    frame : int
        Current frame index

    Returns
    -----------
    edges : list of sets
        Each hydrogen bonded pair (edges of the graph) is returned as a list
    """
    print(frame)

    #Contains number of cells
    cell_num = cell_size ** 3
    
    #Head of cell linked list for all cells
    head_list = np.zeros(cell_num+1)

    #o_list links each atom through the head list in each cell
    o_list = np.zeros(int(water_xyz.shape[0]/3)+1)

    #Fill in the head_list and o_list
    head_list, o_list = make_head(head_list, o_list, water_xyz, box, cell_size, cutoff)

    #Determine hydrogen bonds/edges
    edges = get_frame_edges(water_xyz, box, cell_map, cell_num, head_list, o_list)

    return edges

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

def check_coords(water_xyz, box):
    """
    Parameters 
    -----------
    water_xyz : np.ndarray
        Array containing the positions of all water atoms
    box : float
        Side of the cubic box

    Returns
    -----------
    water_xyz : np.ndarray
        Array with the adjusted positions
    """
    
    #Index of atoms that have position values greater than the box length in the x, y or z directions
    inds = np.where(water_xyz >= box)

    #Move these atom positions back into the box
    water_xyz[inds] -= box

    #Index of atoms that have position values less than 0 in the x, y or z directions
    inds = np.where(water_xyz < 0.0)

    #Move these atom positions back into the box
    water_xyz[inds] += box

    return water_xyz

def get_water_xyz(traj, water_atoms):
    """
    Parameters
    -----------
    traj : MDTraj Trajectory
        system trajectory
    water_atoms : indices of atoms that are in a water molecule

    Returns
    -----------
    water_xyz : np.ndarray
        only contains positions of water atoms
    """
    return traj.xyz[:, water_atoms, :]
    

def get_water_atoms(traj, res_name):
    """
    Parameters
    -----------
    traj : Trajectory object
        MDTraj object
    res_name : str
        Name of the water residues in the PDB File

    Returns
    -----------
    water_residues : list
        list of water molecule residue indices
    water_atoms : list
        list of water atom indices
    """
    water_residues, water_atoms = [], []
    for res in traj.top.residues:
        if res.name == "HOH":
            water_residues.append(res.index)
            for atom in res.atoms:
                water_atoms.append(atom.index)
    return water_residues, water_atoms

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
    water_hist, water_edges, water_p = plt.hist(diams, weights=weights_water, bins=40, alpha=0.5, color='midnightblue',label='Water')
    plt.legend()
    plt.xlabel("Diameter Length")
    plt.ylabel("Probability")

    plt.savefig("diam_h3o_water.png")

def main():
    
    #Determine number of processes you are using (should be 1 per core)
    num_cores = 2

    #Load in MDTraj trajectory
    traj = md.load(sys.argv[1], top=sys.argv[2])

    #Get list of water residues and water atom indices
    water_residues, water_atoms = get_water_atoms(traj, "HOH")

    #List of water residues indexed from 0
    res_index = [i for i in range(len(water_residues))]
    
    frames = len(traj)

    #Get np array (shape = (len(traj), len(water_atoms), 3)) that only contains the water atoms
    water_xyz = get_water_xyz(traj, water_atoms)

    #Make sure coordinates fall between 0 and the box length. Note that this code only works with a cubic box
    water_xyz = check_coords(water_xyz, traj.unitcell_lengths[0,0])

    #Divide the box into cells. The argument specifies the number of cells along each dimension of the box
    cell_map = make_map(10)

    num_frames = len(traj)

    #Call get_network using the joblib Parallel object
    edges = Parallel(n_jobs=num_cores, backend="multiprocessing")(delayed(get_network)(water_xyz[frame], traj.unitcell_lengths[0,0], cell_map, 10, 0.4, frame) for frame in range(num_frames))

    pickle.dump(edges, open("edges.pkl", "wb"))

    ##############################################################################################################################################
    #Comment the above code if you have run the program once so you don't have to get the edges each time you want to try computing a new property
    ##############################################################################################################################################
    edges = pickle.load(open("edges.pkl", "rb"))

    #With the edges obtained from get_network, add to a NetworkX graph object
    graphs = []
    for edge in edges:
        graph = make_graph(edge, res_index)
        graphs.append(graph)
    
    #Uses the formed graphs to compute graph properties
    diams = Parallel(n_jobs=num_cores,backend="multiprocessing")(
            delayed(
                ###Change the "compute_metric" function here if you want to compute a different property
                compute_metric
                )(graphs[frame]) for frame in range(len(graphs)))

    #Change this function depending on which metric you want to plot or analyze
    plt_metric(diams)

if __name__ == "__main__":
    main()
