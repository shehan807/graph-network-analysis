from src.util import *  
import networkx as nx
import pickle
def make_graph(edge, residues):
    graph = nx.Graph()
    for i in residues: graph.add_node(i)
    graph.add_edges_from(edge)
    return graph
def compute_metric(graph):
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

    """
    print(f"processing frame {frame}")
    head_list, linked_list = make_head(xyz, box, num_cells, cutoff, atom_per_res)
    cell_map = make_map(num_cells)
    total_cells = num_cells*num_cells*num_cells 

    edges = []
    # Get list of "RESNAME.ATOMNAME" from traj_filtered (should be length atom_per_residue*len(residues))
    # filtered_atoms = [f"{atom.residue.name}.{atom.name}" for atom in traj_filtered.topology.atoms]
    filtered_atoms = [f"{atom.name}" for atom in traj_filtered.topology.atoms]
    print(f"filtered atom list length: {len(filtered_atoms)}")
    #print(filtered_atoms)
    print(len(residues),len(atoms))
    for i in range(total_cells):
        i_mol = head_list[i] # access initial molecule from head list 
        # print(i_mol)
        #print(f"found i_mol={i_mol} from head_list[i].")
        while i_mol > 0: # follow trail of link-list references until element is zero
            j_mol = linked_list[int(i_mol)]
            while j_mol > 0:
                i_mol_atom_ind = int(i_mol*atom_per_res) 
                j_mol_atom_ind = int(j_mol*atom_per_res) 
                # obtain slice of filtered_atoms for just i_mol and j_mol
                i_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                j_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                #print(i_mol, j_mol)
                #print(i_mol_atom_ind, j_mol_atom_ind)
                #print(i_mol_atoms, j_mol_atoms)
                if check_criteria(criteria, filtered_atoms, i_mol_atom_ind, i_mol_atoms, j_mol_atom_ind, j_mol_atoms, xyz, traj_filtered): edges.append((int(i_mol), int(j_mol)))
                j_mol = linked_list[int(j_mol)]
        
        
            # loop through neighboring cells 
            jcell0 = i*13 # double check that this is correct, should be...
            #print(f"jcell0 = {jcell0}; i = {i}")
            for neighbor in range(13):
                #print(f"\n %%% Neighbor {neighbor} %%% \n")
                # print(f"jcell0 = {jcell0}; i = {i}; cell_map ind = {jcell0 + neighbor}")
                jcell = int(cell_map[jcell0 + neighbor])
                # head molecule index of jcell 
                j_mol = head_list[jcell]
                
                while j_mol > 0: 
                    i_mol_atom_ind = int(i_mol*atom_per_res) 
                    j_mol_atom_ind = int(j_mol*atom_per_res) 
                    # obtain slice of filtered_atoms for just i_mol and j_mol
                    i_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                    j_mol_atoms = filtered_atoms[i_mol_atom_ind:i_mol_atom_ind+atom_per_res]
                    #print(i_mol, j_mol)
                    #print(i_mol_atom_ind, j_mol_atom_ind)
                    #print(i_mol_atoms, j_mol_atoms)
                    if check_criteria(criteria, filtered_atoms, i_mol_atom_ind, i_mol_atoms, j_mol_atom_ind, j_mol_atoms, xyz, traj_filtered): edges.append((int(i_mol), int(j_mol)))

                    j_mol = linked_list[int(j_mol)]
            i_mol = linked_list[int(i_mol)]

    return edges
