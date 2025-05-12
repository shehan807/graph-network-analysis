#!/bin/bash
# EXDIR="" # running from example dir
EXDIR="examples/PTFA" # running from main dir
#python3 main.py\
#	--config_file "$EXDIR"/config.yaml\
#	--dcd_file "$EXDIR"/inputs/md_npt-4_2p5.dcd\
#	--pdb_file "$EXDIR"/inputs/npt_final-4_2p5.pdb

python3 main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/md_npt-4_5p0.dcd\
	--pdb_file "$EXDIR"/inputs/npt_final-4_5p0.pdb\
	--output_directory "$EXDIR"/outputs/0.05	

#python3 main.py\
#	--config_file "$EXDIR"/config.yaml\
#	--dcd_file "$EXDIR"/inputs/md_npt-4_7p5.dcd\
#	--pdb_file "$EXDIR"/inputs/npt_final-4_7p5.pdb\
#	--output_directory "$EXDIR"/outputs/0.075	

python3 main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/md_npt-4_10p0.dcd\
	--pdb_file "$EXDIR"/inputs/npt_final-4_10p0.pdb\
	--output_directory "$EXDIR"/outputs/0.1	
