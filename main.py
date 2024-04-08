import os
import argparse 
import yaml
from pathlib import Path
from src import visualization, util
from src.network_analysis import * 
from joblib import Parallel, delayed 
import time
import pickle

def main():
    parser = argparse.ArgumentParser(description="Graph Network Analysis for Molecular Dynamics Trajectories")
    parser.add_argument("-i", '--config_file', type=str, help='Path to the YAML configuration file', required=True)
    parser.add_argument("-d", '--dcd_file',  type=str, help='Path to the input .dcd MD trajectory file', required=True)
    parser.add_argument("-p", '--pdb_file',  type=str, help='Path to the input .pdb MD topology file', required=True)
    parser.add_argument("-o", '--output_directory', type=str, help='Path to the output directory', required=True)
    
    args = parser.parse_args()

    with open(args.config_file, 'r') as config_file:
        config = yaml.safe_load(config_file)
    
    # create path objects 
    inp = Path(args.config_file)
    dcd = Path(args.dcd_file)
    pdb = Path(args.pdb_file)
    out = Path(args.output_directory)
    if not os.path.exists(out): os.makedirs(out)
    
    # obtain config.yaml variables
    num_cores = config["num_cores"]
    num_cells = config["num_cells"]
    cutoff = config["cutoff"]
    residue_name = config["residue_name"]
    criteria = config["criteria"]
    check_pkl = config["check_pkl"]
    

    # obtain network prerequisites
    xyz, traj, traj_filtered, residues, atoms, num_frames, box, atom_per_res = util.initialize(dcd, pdb, residue_name, num_cells, cutoff)
    res_index = [i for i in range(len(residues))]
    #num_frames = 1 
    # obtain edges 
    #if os.environ["MULTINODE"] == True:
    #    pass
    #else:
    edges = Parallel(n_jobs=num_cores,backend="multiprocessing")(delayed(get_network)(xyz[frame], box, num_cells, cutoff, frame, atom_per_res, residues, atoms, traj, traj_filtered, criteria) for frame in range(num_frames))
    pickle.dump(edges, open(os.path.join(out,"edges.pkl"), "wb"))
    
    with open("edges.txt", "w") as etxt:
        for edge in edges[0]:
            etxt.writelines(str(edge)+"\n")

    # save or load contents of edges for later use
    # edges = pickle.load(open(os.path.join(out,"edges.pkl"), "rb"))

    # use edges obtained from get_network, add to a NetworkX graph object
    graphs = []
    for edge in edges:
        graph = make_graph(edge, res_index)
        graphs.append(graph)
     
    print(f"graph len: {len(graphs)}")
    # use the formed graphs to compute graph properties
    diams = Parallel(n_jobs=1,backend="multiprocessing")(
            delayed(
                ###Change the "compute_metric" function here if you want to compute a different property
                compute_metric
                )(graphs[frame]) for frame in range(len(graphs)))
   
    
    with open("diam.txt", "w") as dtxt:
        for diam in diams: 
            #for d in diam: 
            dtxt.writelines(str(diam)+"\n")
    # change this function depending on which metric you want to plot or analyze
    visualization.plt_metric(diams, out)

if __name__ == "__main__":
    main()
