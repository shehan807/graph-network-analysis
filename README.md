# Graph Network Code for Analyzing Molecular Dynamics Trajectories
Code for calculating inter-species cluster network diameters for various MD simulations. 
# Contents
* `src/` core code
    * `network_analysis.py` functions for creating graphs and calculating graph metrics
    * `visualization.py` plotting diameters of graph network
    * `util.py` functions for computing distances, angles, and creating cell maps
    	* `initialize`: 
    	* `get_atoms`: 
    	* `make_head`, `set_cell`, `make_map`: 
* `examples/` 
    * `water/` See  Int. J. Mol. Sci. 2020, 21(2), 403; https://doi.org/10.3390/ijms21020403 
    * `N1888/` 
# Example: Water cluster calculations
## Installation 
``` git clone ...```
## Setting up confing.yaml
```
config file 
```
## Running shell script

# Copyright
Copyright (c) 2024, McDaniel Group at Georgia Institute of Technology and and Authors:
Authors:
- [Shehan Parmar](https://shehan807.github.io/) (sparmar32@gatech.edu
- John P. Stoppleman (jstoppelman3@gatech.edu)
