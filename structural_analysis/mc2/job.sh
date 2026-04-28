#!/bin/bash
set -euo pipefail

DIRS=("mc1" "mc2" "mc3")

for dir in "${DIRS[@]}"; do
    echo "Processing $dir"

    cd "$dir" || exit

    for e in 2 4; do
        echo "  Running with e=$e"

        # Replace <e> with value and save as cn.m
        sed "s/<e>/$e/g" cn0.m > cn.m

        # Run MATLAB in non-interactive mode
	matlab -nodisplay -r "run partial_rdf.m; exit()" >> check.dat

        # Optional: save outputs separately
        mv check.dat "check_${e}.txt" 2>/dev/null || true
    done

    cd ..
done
