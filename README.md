# Graph Network Code for Analyzing Molecular Dynamics Trajectories
Code for calculating inter-species cluster network diameters for various MD simulations. 
# Contents
* `src/` core code
    * `network_analysis.py` functions for creating graphs and calculating graph metrics
    	* `make_graph`: create NetworkX Graph
    	* `compute_metric`: calculate NetworkX Graph [diameter](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.distance_measures.diameter.html) 
    	* `get_network`: apply [cell-linked list algorithm](https://github.com/glennklockwood/allen-tildesley/blob/master/link-cell.f) to generate edge list 
    * `visualization.py` plotting diameters of graph network
    * `util.py` functions for computing distances, angles, and creating cell maps
* `examples/` 
    * `water/` See  [Int. J. Mol. Sci. 2020, 21(2), 403](https://doi.org/10.3390/ijms21020403) 
    * `N1888/` See [J. Phys. Chem. B 2024.](https://doi.org/10.1021/acs.jpcb.4c06255)
      
# Example: Water cluster calculations
Reproducing the water cluster analysis results requires 3 steps. Following the `examples/water' directory: 

## 1. Installation 
``` 
git clone https://github.com/shehan807/graph-network-analysis.git
```

## 2. Set up confing.yaml
The `config.yaml` contains all tunable parameters (number of cells, cutoff distances, edge criteria definitions, etc) to obtain a cluster of like-species from an MD trajectory. The example below specifies a water molecule residue, `HOH` along with a distance criteria of 0.36 nanometers to generate two oxygen atoms, `O-O`, as nodes connected by an edge. 

```yaml
num_cores: 1 # Number of cores/processes for joblib.Parallel 
num_cells: 10 # Number of cells to divide along each dimension of MD box 
cutoff: 0.4 # cutoff distance, nm
residue_name: "HOH" # species for inter-species cluster analysis
criteria: # list of criteria for network edge
  - name: "O-O" # atom symbols 
    distance: 0.36 # nm 
    min_true: 1
    angle: None # deg
      #  - name: "H1-O"
      #    distance: 0.245
      #    angle: 30.0
      #  - name: "H2-O"
      #    distance: 0.245
      #    angle: 30.0
check_pkl: True # check if edges.pkl exists in code
```
For automating the process of creating numerous `config.yaml` files, see `examples/N1888/create_config.sh`.

## 3. Run Shell Script
Create an `inputs` directory to place a `<topology>.pdb` and `<trajectory>.dcd` file and run

```
./run_analysis
```

The `/outputs` directory will contain a pickle file for the edges (to avoid having to rerun NetworkX edge generation in the future) along with the generated probability ditribution of diameters, $P(d)$, that should look something like:

![alt text](https://github.com/shehan807/graph-network-analysis/blob/main/examples/N1888/outputs/diam.png)

# Copyright
Copyright (c) 2024, McDaniel Group at Georgia Institute of Technology and and Authors:
Authors:
- [Shehan Parmar](https://shehan807.github.io/) (sparmar32@gatech.edu)
- John P. Stoppleman (jstoppelman3@gatech.edu)
