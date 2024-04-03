import argparse 
import yaml
from pathlib import Path
from src import network_analysis, visualization, util

def main():
    parser = argparse.ArgumentParser(description="Graph Network Analysis for Molecular Dynamics Trajectories")
    parser.add_argument("-i", '--config_file', type=str, help='Path to the YAML configuration file', required=True)
    parser.add_argument("-d", '--dcd_file',  type=str, help='Path to the input .dcd MD trajectory file', required=True)
    parser.add_argument("-p", '--pdb_file',  type=str, help='Path to the input .pdb MD topology file', required=True)
    parser.add_argument("-o", '--output_directory', type=str, help='Path to the output directory', required=True)
    
    args = parser.parse_args()

    with open(args.config_file, 'r') as config_file:
        config = yaml.safe_load(config_file)
    
    inp = Path(args.config_file)
    dcd = Path(args.dcd_file)
    pdb = Path(args.pdb_file)
    out = Path(args.output_directory)
    num_cores = config["num_cores"]
    num_cells = config["num_cells"]
    cutoff = config["cutoff"]
    residues = config["residue_names"]
    criteria = config["criteria"]
    check_pkl = config["check_pkl"]




if __name__ == "__main__":
    main()
