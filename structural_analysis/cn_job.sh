#!/bin/bash
set -euo pipefail

DIRS=("mc1" "mc2" "mc3" "mc4")

for dir in "${DIRS[@]}"; do
    echo "Processing $dir"

    cd "$dir" || exit

    for e in 2 4; do
        echo "  Running with e=$e"

        # Assign r based on e
        if [[ "$e" -eq 2 ]]; then
            r=2.35
        else
            r=2.55
        fi

        # Replace both <e> and <r> in one go
        sed -e "s/<e>/$e/g" -e "s/<r>/$r/g" cn_0.m > cn.m

        # Run MATLAB in non-interactive mode
        matlab -nodisplay -nosplash -r "run('cn.m'); exit;" >> temp.txt 2>&1

        # Save output separately
        mv temp.txt "cn_${e}.txt" 2>/dev/null || true
    done

    cd ..
done
