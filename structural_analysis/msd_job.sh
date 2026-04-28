#!/bin/bash
set -euo pipefail

DIRS=("mc1" "mc2" "mc3" "mc4")

for dir in "${DIRS[@]}"; do
    echo "Processing $dir"

    cd "$dir" || exit

    echo "  Running with e=$dir"


    # Run MATLAB in non-interactive mode
    matlab -nodisplay -r "run self_diffusivity.m; exit()" >> temp.txt

    # save outputs separately
    mv temp.txt "msd.txt" 2>/dev/null || true
    
    cd ..
done
