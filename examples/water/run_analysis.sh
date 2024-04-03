#!/bin/bash
# EXDIR="" # running from example dir
EXDIR="examples/water" # running from main dir
python3 main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/out.dcd\
	--pdb_file "$EXDIR"/inputs/equilibrated.pdb\
	--output_directory "$EXDIR"/outputs	
