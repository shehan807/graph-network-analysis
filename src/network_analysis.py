from src.util import *  
import networkx as nx
import pickle
def make_graph():
    pass
def get_network(xyz, box, num_cells, cutoff, frame, atom_per_res, residues, atoms, traj, criteria):
    """
    Parameters
    ----------
    xyz : np.ndarray
        Array containing positions of selected residue 

    """
    print(f"processing frame {frame}")
    head_list, linked_list = make_head(xyz, box, num_cells, cutoff, atom_per_res)
    cell_map = make_map(num_cells)
   
    edges = []
    for i in range(num_cells):
        i_mol = head_list[i] # access initial molecule from head list 
        print(f"found i_mol={i_mol} from head_list[i].")
        while i_mol > 0: # follow trail of link-list references until element is zero
            j_mol = linked_list[int(i_mol)]
            while j_mol > 0:
                j_mol = linked_list[int(j_mol)]
                # get list of "RES.ATOM" corresponding to traj_filtered (should be length atom_per_residue*len(residues))
                # get slice of RESATOM list corresponding to i_mol and j_mol 
                # criteria_bool = np.zeros(len(criteria))
                # create temp edge list for each criteria
                # for criterion in criteria:
                    # create pair list of evaluated pairs, i.e., ["ATOM-ATOM",...]
                    # create unique atom list of atoms evaluated, i.e. ["ATOM", "ATOM", ...]
                    # create list of atom indices from sliced RESATOM list that are in unique atom list 
                    # for i_mol_matched_atoms: 
                        # for j_mol_match_atoms: 
                            # evaluate distance criterion
                            # evaluate angle criterion
                            # update criteria_bool 
                # if all criteria are true: edges.append(edge)


               

    return edges
def compute_metric(graph):
    """
    """
    pass
