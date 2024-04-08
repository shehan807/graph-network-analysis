#!/bin/bash
values=(C00 C01 C03 C04 C05 C06 C07 C09 C0A C0B C0C C0D C0E C0F C0G C0I C0J C0K C0M C0N C0O C0P C0Q)
outfile="combinations.txt"
rm "$outfile"
# Loop over each element in the array
for i in "${values[@]}"; do
    for j in "${values[@]}"; do
        # Print the configuration lines for each pair
        echo "  - name: \"$i-$j\" # atom symbols" >> "$outfile"
        echo "    distance: 0.70 # nm" >> "$outfile"
    done
done
