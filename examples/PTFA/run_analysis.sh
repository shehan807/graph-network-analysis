#!/bin/bash
# EXDIR="" # running from example dir
EXDIR="." # examples/PTFA" # running from main dir
### Temp = 300K 
python3 ../../main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/300K/md_npt-4_2p5.dcd\
	--pdb_file "$EXDIR"/inputs/300K/npt_final-4_2p5.pdb\
	--output_directory "$EXDIR"/outputs/300K/0.025	

python3 ../../main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/300K/md_npt-4_5p0.dcd\
	--pdb_file "$EXDIR"/inputs/300K/npt_final-4_5p0.pdb\
	--output_directory "$EXDIR"/outputs/300K/0.05	

python3 ../../main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/300K/md_npt-4_7p5.dcd\
	--pdb_file "$EXDIR"/inputs/300K/npt_final-4_7p5.pdb\
	--output_directory "$EXDIR"/outputs/300K/0.075	

python3 ../../main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/300K/md_npt-4_10p0.dcd\
	--pdb_file "$EXDIR"/inputs/300K/npt_final-4_10p0.pdb\
	--output_directory "$EXDIR"/outputs/300K/0.1
### Temp = 373K	
python3 ../../main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/373K/md_npt-4_2p5.dcd\
	--pdb_file "$EXDIR"/inputs/373K/npt_final-4_2p5.pdb\
	--output_directory "$EXDIR"/outputs/373K/0.025	

python3 ../../main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/373K/md_npt-4_5p0.dcd\
	--pdb_file "$EXDIR"/inputs/373K/npt_final-4_5p0.pdb\
	--output_directory "$EXDIR"/outputs/373K/0.05	

python3 ../../main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/373K/md_npt-4_7p5.dcd\
	--pdb_file "$EXDIR"/inputs/373K/npt_final-4_7p5.pdb\
	--output_directory "$EXDIR"/outputs/373K/0.075	

python3 ../../main.py\
	--config_file "$EXDIR"/config.yaml\
	--dcd_file "$EXDIR"/inputs/373K/md_npt-4_10p0.dcd\
	--pdb_file "$EXDIR"/inputs/373K/npt_final-4_10p0.pdb\
	--output_directory "$EXDIR"/outputs/373K/0.1	
