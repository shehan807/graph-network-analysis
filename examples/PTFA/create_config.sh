#!/bin/bash
values=(C00 C01 C03 C04 C05 C06 C07 C09 C0A C0B C0C C0D C0E C0F C0G C0I C0J C0K C0M C0N C0O C0P C0Q)
outfile="combinations.txt"
rm "$outfile"

echo "num_cores: 1 # Number of cores/processes for joblib.Parallel" >> "$outfile" 
echo "num_cells: 10 # Number of cells to divide along each dimension of MD box" >> "$outfile"
echo "cutoff: 0.4 # cutoff distance, nm" >> "$outfile"
echo 'residue_name: "OCT" # species for inter-species cluster analysis' >> "$outfile"
echo "criteria: # list of criteria for network edge" >> "$outfile"

# Loop over each element in the array
for i in "${values[@]}"; do
    for j in "${values[@]}"; do
        # Print the configuration lines for each pair
        echo "  - name: \"$i-$j\" # atom symbols" >> "$outfile"
        echo "    distance: 0.45 # nm" >> "$outfile"
        echo "    min_true: 2 " >> "$outfile"
    done
done

echo "check_pkl: True # check if edges.pkl exists in code" >> "$outfile"
