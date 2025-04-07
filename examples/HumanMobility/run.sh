
#!/bin/bash

# Remove previously generated data and images to regenerate the new ones
rm humans.hdf5
rm river.hdf5
rm data/*
#rm images/*

# Generate the initial conditions if they are not present.
if [ ! -e humans.hdf5 ]
then
    echo "Generating initial conditions for the human mobility box example..."
    python3 makeIC.py -t gas -f humans.hdf5  # -t particles
fi

# Generate acceleration field for river
if [ ! -e river.hdf5 ]
then
    echo "Generating acceleration field for the river..."
    python3 makeRiver.py -f river.hdf5
fi

# Run SWIFT
# gdb --args swift -g --threads=4 -n 10000 humanMobility.yml # -A -s
#time swift --hm-river --hm-randomwalk --threads=8 -n 10000 humanMobility.yml # -A -s -g -G 
time mpirun -n 2 swift_mpi -A -s -g -G --hm-river --hm-randomwalk --threads=4 -n 50000 humanMobility.yml #
