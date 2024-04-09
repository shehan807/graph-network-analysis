from src.util import *  
import networkx as nx
import pickle
import logging 

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
    diam = []
    clusters = [c for c in sorted(nx.connected_components(graph),key=len,reverse=True)]
    clusters_lengths = [len(c) for c in sorted(nx.connected_components(graph),key=len,reverse=True)]
    logging.info("Found %d clusters (covering %d total nodes/molecules)",len(clusters), sum(clusters_lengths))
    logging.info("cluster_lengths list = %s",clusters_lengths)

    for cluster in nx.connected_components(graph):
        diam.append(nx.diameter(graph.subgraph(cluster)))
    return diam

def get_network(xyz, box, num_cells, cutoff, frame, atom_per_res, residues, atoms, traj, out, criteria):
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
    
    # because traj_filtered is read_only, joblib.Parallel can only work this way
    trj_fil_file = os.path.join(out, "trj_filtered.pkl")
    if os.path.isfile(trj_fil_file):
        with open(trj_fil_file, "rb") as file:
            traj_filtered = pickle.load(file)
    else: 
        traj_filtered = traj.atom_slice(atoms)
        with open(trj_fil_file, "wb") as file: 
            pickle.dump(traj_filtered, file)

    head_list, linked_list = make_head(xyz, box, num_cells, cutoff, atom_per_res)
    cell_map = make_map(num_cells)
    total_cells = num_cells*num_cells*num_cells 

    edges = []
    # Get list of "ATOMNAMES" from traj_filtered (length atom_per_residue*len(residues))
    filtered_atoms = [f"{atom.name}" for atom in traj_filtered.topology.atoms]
    for i in range(total_cells):
        i_mol = head_list[i] # access initial molecule from head list 
        #Some cells won't have a molecule, so only need to loop through the cell if i_mol isn't 0
        i_mol_cond = True
        while i_mol > 0: # follow trail of link-list references until element is zero
            #if i_mol == 0: i_mol_cond = False
            
            #Get the linked j_atom from linked_list
            j_mol = linked_list[int(i_mol)]
            
            j_mol_cond = True
            while j_mol > 0:
                
                # set end condition to be when j_mol reach end of linked_list trail
                #if j_mol == 0: j_mol_cond = False
                #if (i_mol == 0) and (j_mol == 0): continue
                
                # index of atom from i_mol and j_mol 
                i_mol_atom_ind = int(i_mol*atom_per_res) 
                j_mol_atom_ind = int(j_mol*atom_per_res) 
                
                # obtain slice of filtered_atoms for just i_mol and j_mol
                i_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                j_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                
                # check all criteria
                if check_criteria(criteria, filtered_atoms, i_mol, i_mol_atom_ind, i_mol_atoms, j_mol, j_mol_atom_ind, j_mol_atoms, xyz, traj_filtered): edges.append((int(i_mol), int(j_mol)))
                j_mol = linked_list[int(j_mol)]
        
            # loop through neighboring cells 
            jcell0 = i*13 # double check that this is correct, should be...
            for neighbor in range(13):
                jcell = int(cell_map[jcell0 + neighbor])
                
                # head molecule index of jcell 
                j_mol = head_list[jcell]
                
                j_mol_cond = True
                while j_mol > 0: 
                    
                    # set end condition to be when j_mol reach end of linked_list trail
                    #if j_mol == 0: j_mol_cond = False
                    #if (i_mol == 0) and (j_mol == 0): continue
                    
                    # index of atom from i_mol and j_mol
                    i_mol_atom_ind = int(i_mol*atom_per_res) 
                    j_mol_atom_ind = int(j_mol*atom_per_res) 
                    
                    # obtain slice of filtered_atoms for just i_mol and j_mol
                    i_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                    j_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                    # check all criteria again
                    if check_criteria(criteria, filtered_atoms, i_mol, i_mol_atom_ind, i_mol_atoms, j_mol, j_mol_atom_ind, j_mol_atoms, xyz, traj_filtered): edges.append((int(i_mol), int(j_mol)))

                    #See if there's an atom linked to j_mol
                    j_mol = linked_list[int(j_mol)]
            
            #See if there's a linked molecule to i_mol
            i_mol = linked_list[int(i_mol)]

    return edges
