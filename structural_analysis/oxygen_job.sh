#!/bin/bash
set -euo pipefail

DIRS=("mc1" "mc2" "mc3" "mc4")

for dir in "${DIRS[@]}"; do
    echo "Processing $dir"

    cd "$dir" || exit

    echo "  Running with e=$dir"


    # Run MATLAB in non-interactive mode
    matlab -nodisplay -r "run oxygen_type.m; exit()" >> oxygen_type.txt

    cd ..
done
