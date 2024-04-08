#!/bin/bash
# EXDIR="" # running from example dir
EXDIR="examples/N1888" # running from main dir
python3 main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/npt_final-28.dcd\
	--pdb_file "$EXDIR"/inputs/npt_final-28.pdb\
	--output_directory "$EXDIR"/outputs	
