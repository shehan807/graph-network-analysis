import os
import argparse 
import yaml
from pathlib import Path
from src import visualization, util
from src.network_analysis import * 
from joblib import Parallel, delayed 
import time
import pickle
import logging

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
    
    # configure logging in python 
    logging.basicConfig(filename=os.path.join(out,'summary.out'),level=logging.DEBUG, format='%(levelname)s : %(message)s', filemode='w')

    # obtain config.yaml variables
    num_cores = config["num_cores"]
    num_cells = config["num_cells"]
    cutoff = config["cutoff"]
    residue_name = config["residue_name"]
    criteria = config["criteria"]
    check_pkl = config["check_pkl"]

    # obtain network prerequisites
    xyz, traj, traj_filtered, residues, atoms, num_frames, box, atom_per_res = util.initialize(dcd, pdb, out, residue_name, num_cells, cutoff)
    res_index = [i for i in range(len(residues))]
    logging.info("%d %s residues found.", len(res_index), residue_name)

    # num_frames = 50 
    # obtain edges 
    #if os.environ["MULTINODE"] == True:
    #    pass
    #else:
    start = time.time()
    edges = Parallel(n_jobs=num_cores,backend="multiprocessing")(delayed(get_network)(xyz[frame], box, num_cells, cutoff, frame, atom_per_res, residues, atoms, traj, out, criteria) for frame in range(num_frames))
    end = time.time()
    total_time = (end - start) / 60 
    logging.info("%d frames took %.2f min", num_frames, total_time)
    logging.info("edges:\n%s", edges)

    pickle.dump(edges, open(os.path.join(out,"edges.pkl"), "wb"))
    
    # save or load contents of edges for later use
    # edges = pickle.load(open(os.path.join(out,"edges.pkl"), "rb"))

    # use edges obtained from get_network, add to a NetworkX graph object
    graphs = []
    for edge in edges:
        graph = make_graph(edge, res_index)
        graphs.append(graph)
     
    # use the formed graphs to compute graph properties
    diams = Parallel(n_jobs=num_cores,backend="multiprocessing")(
            delayed(
                ###Change the "compute_metric" function here if you want to compute a different property
                compute_metric
                )(graphs[frame]) for frame in range(len(graphs)))
   
    # change this function depending on which metric you want to plot or analyze
    visualization.plt_metric(diams, out)

if __name__ == "__main__":
    main()
