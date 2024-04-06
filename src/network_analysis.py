from src.util import *  
import networkx as nx
import pickle

def make_graph(edge, res_index):
    """
    Parameters
    -----------
    edge : list
        contains pairs of connected molecule nodes
    res_index : list
        index of each molecule

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
        List of each graph cluster diameter
    """
    # compute diameter of the cluster here
    diam_water = []
    for cl in nx.connected_components(graph):
        diam_water.append(nx.diameter(graph.subgraph(cl)))
    return diam_water

def get_network(xyz, box, num_cells, cutoff, frame, atom_per_res, residues, atoms, traj, traj_filtered, criteria):
    """
    Parameters
    ----------
    xyz : np.ndarray
        Array containing positions of selected residue
    box : float 
        trajectory object box dimensions
    num_cells : float 
        number of cells in one dimension 
    cutoff : float
        The box dimension divided by the cell_size must be larger than this amount. Don't adjust this as it should probably be larger than the O-O hydrogen bond distance
    frame : int
        Current frame index
    atom_per_res : int
        number of atoms per selected molecule (e.g., H2O = 3)
    residues : list
        list of selected residues
    atoms : list 
        list of atoms corresponding to residue list
    traj : mdtraj.core.trajectory.Trajectory 
        MDTraj trajectory object
    traj_filtered : mdtraj.core.trajectory.Trajectory 
        MDTraj trajectory object of only relevant residues
    criteria : dict
        dictionary of criteria read in from config.yaml
    """
    print(f"processing frame {frame}")
    head_list, linked_list = make_head(xyz, box, num_cells, cutoff, atom_per_res)
    cell_map = make_map(num_cells)
    total_cells = num_cells*num_cells*num_cells 

    edges = []
    # Get list of "ATOMNAMES" from traj_filtered (length atom_per_residue*len(residues))
    filtered_atoms = [f"{atom.name}" for atom in traj_filtered.topology.atoms]
    for i in range(total_cells):
        i_mol = head_list[i] # access initial molecule from head list 
        #Some cells won't have a molecule, so only need to loop through the cell if i_mol isn't 0
        while i_mol > 0: # follow trail of link-list references until element is zero
            #Get the linked j_atom from linked_list
            j_mol = linked_list[int(i_mol)]
            while j_mol > 0:
                
                # index of atom from i_mol and j_mol 
                i_mol_atom_ind = int(i_mol*atom_per_res) 
                j_mol_atom_ind = int(j_mol*atom_per_res) 
                
                # obtain slice of filtered_atoms for just i_mol and j_mol
                i_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                j_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                
                # check all criteria
                if check_criteria(criteria, filtered_atoms, i_mol_atom_ind, i_mol_atoms, j_mol_atom_ind, j_mol_atoms, xyz, traj_filtered): edges.append((int(i_mol), int(j_mol)))
                j_mol = linked_list[int(j_mol)]
        
            # loop through neighboring cells 
            jcell0 = i*13 # double check that this is correct, should be...
            for neighbor in range(13):
                jcell = int(cell_map[jcell0 + neighbor])
                
                # head molecule index of jcell 
                j_mol = head_list[jcell]
                
                while j_mol > 0: 
                    
                    # index of atom from i_mol and j_mol
                    i_mol_atom_ind = int(i_mol*atom_per_res) 
                    j_mol_atom_ind = int(j_mol*atom_per_res) 
                    
                    # obtain slice of filtered_atoms for just i_mol and j_mol
                    i_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                    j_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                    # check all criteria again
                    if check_criteria(criteria, filtered_atoms, i_mol_atom_ind, i_mol_atoms, j_mol_atom_ind, j_mol_atoms, xyz, traj_filtered): edges.append((int(i_mol), int(j_mol)))

                    #See if there's an atom linked to j_mol
                    j_mol = linked_list[int(j_mol)]
            
            #See if there's a linked molecule to i_mol
            i_mol = linked_list[int(i_mol)]

    return edges
