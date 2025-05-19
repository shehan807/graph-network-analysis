#!/bin/bash
#SBATCH -A hive-jmcdaniel43
#SBATCH --job-name=water-cluster
#SBATCH --time=24:00:00
#SBATCH -N1 -n16
#SBATCH --partition=hive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=parmar@gatech.edu
#SBATCH --output=energies.%x-%j
#SBATCH --error=error.%x-%j

hostname
module load anaconda3

# fill this in with your OpenMM conda environment
source activate cluster-analysis

bash run_analysis.sh &>> output.txt

end=`date +%s`
echo "Done   " `date`
runtime=$((end-start))
hours=$(printf %02d $((runtime / 3600)))
minutes=$(printf %02d $(( (runtime % 3600) / 60 )))
seconds=$(printf %02d $(( (runtime % 3600) % 60 )))
echo "Elapsed time: $hours:$minutes:$seconds (hh:mm:ss)"
echo "Elapsed time: $runtime seconds"
exit $ret
