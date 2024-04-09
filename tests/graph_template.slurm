#!/bin/bash
#SBATCH --account=gts-jmcdaniel43-chemx
#SBATCH --job-name=${VAR1}_${VAR2}_cluster
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=01:00:00
#SBATCH --partition=inferno
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=parmar@gatech.edu
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out

hostname

module purge
module load anaconda3

source activate cluster-analysis

distance_value=${VAR1} # First input: distance value
min_true_value=${VAR2}  # Second input: min_true value

./create_config.sh $distance_value $min_true_value $SLURM_CPUS_ON_NODE
cp combinations.txt config_"$distance_value"_"$min_true_value".yaml

EXDIR="/storage/home/hcoda1/4/sparmar32/p-jmcdaniel43-0/N1888/pol/temp-study/2x2x2/300K/ANALYSIS/graph-cluster/graph-network-analysis/tests/N1888" # running from main dir
CODEDIR="/storage/home/hcoda1/4/sparmar32/p-jmcdaniel43-0/N1888/pol/temp-study/2x2x2/300K/ANALYSIS/graph-cluster/graph-network-analysis" # running from main dir
DATADIR="/storage/home/hcoda1/4/sparmar32/p-jmcdaniel43-0/N1888/pol/temp-study/2x2x2/300K/simulation_output" # running from main dir

python3 "$CODEDIR"/main.py\
	--config_file "$EXDIR"/config_"$distance_value"_"$min_true_value".yaml\
	--dcd_file "$DATADIR"/md_npt-28.dcd\
	--pdb_file "$DATADIR"/npt_final-28.pdb\
	--output_directory "$EXDIR"/outputs_"$distance_value"_"$min_true_value"	

end=`date +%s`
echo "Done   " `date`
runtime=$((end-start))
hours=$(printf %02d $((runtime / 3600)))
minutes=$(printf %02d $(( (runtime % 3600) / 60 )))
seconds=$(printf %02d $(( (runtime % 3600) % 60 )))
echo "Elapsed time: $hours:$minutes:$seconds (hh:mm:ss)"
echo "Elapsed time: $runtime seconds"
exit $ret


































































