#!/bin/bash

# Define ranges for your variables
var1_values=(0.38 0.385 0.39 0.395 0.400 0.405 0.41 0.4105 0.412 0.413 0.414 0.415 0.416 0.417 0.418 0.419 0.42 0.430 0.44 0.45 0.455 0.46 0.47 0.48) # distance criterion
var2_values=(1 2 3 4 5 6 ) # min_true_vals

# Loop over each combination of variables
count=0
for var1 in "${var1_values[@]}"; do
    for var2 in "${var2_values[@]}"; do
        # Replace placeholders in the template with actual values
        sed "s/\${VAR1}/$var1/g; s/\${VAR2}/$var2/g" graph_template.slurm > graph_${count}.temp
        ((count++))
        # Submit job
        # sbatch graph_"$count".temp
        
        # Optional: wait a bit to avoid overwhelming the scheduler
        sleep 1
    done
done

