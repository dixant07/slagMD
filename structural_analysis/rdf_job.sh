#!/bin/bash
set -euo pipefail

DIRS=("mc1" "mc2" "mc3" "mc4")

for dir in "${DIRS[@]}"; do
    echo "Processing $dir"

    cd "$dir" || exit

    for e in 2 4; do
        echo "  Running with e=$e"

        # Replace <e> with value and save as cn.m
        sed "s/<e>/$e/g" partial_rdf_0.m > partial_rdf.m

        # Run MATLAB in non-interactive mode
	matlab -nodisplay -r "run partial_rdf.m; exit()" >> temp.txt

        # Optional: save outputs separately
        mv temp.txt "rdf_${e}.txt" 2>/dev/null || true
    done

    cd ..
done
