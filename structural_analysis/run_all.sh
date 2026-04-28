#!/bin/bash

# List of folders
FOLDERS=("mc1" "mc2" "mc3" "mc4")  # add as many as you have

# LAMMPS executable with MPI
LAMMPS_EXEC="mpirun -np 12 ../lammps/build/lmp"  # adjust -np for number of cores

# Name of the input script in each folder (assume same name)
INPUT_SCRIPT="slag.in"  # change to your input script name

# Loop over folders
for DIR in "${FOLDERS[@]}"; do
    echo "Entering folder $DIR ..."
    cd "$DIR" || { echo "Failed to enter $DIR"; exit 1; }

    echo "Running LAMMPS in $DIR ..."
    $LAMMPS_EXEC -in "$INPUT_SCRIPT"

    echo "Finished $DIR"
    cd ..  # go back to parent folder
done

echo "All folders completed."
