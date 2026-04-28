#!/bin/bash
set -euo pipefail

DIRS=("1863" "1893" "1923")

for temp in ${DIRS[@]}; do
    echo "Processing $temp"

    cd "$temp" || exit

    echo "  Running with e=$temp"

    sed "s/<t>/$temp/g" self_diffusivity_0.m > self_diffusivity.m
    # Run MATLAB in non-interactive mode
    matlab -nodisplay -r "run self_diffusivity.m; exit()" > temp.txt

    # save outputs separately
    mv temp.txt "msd.txt" 2>/dev/null || true
    
    cd ..
done
