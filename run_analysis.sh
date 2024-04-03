#!/bin/bash
python3 main.py\
	--config_file examples/water/config.yaml\
	--dcd_file examples/water/inputs/out.dcd\
	--pdb_file examples/water/equilibrated.pdb\
	--output_directory examples/water/outputs	
