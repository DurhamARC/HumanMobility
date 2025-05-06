#!/bin/bash

# Remove previously generated data and images to regenerate the new ones
rm humans.hdf5
rm river.hdf5
rm data/*
rm images/*

# Generate the initial conditions if they are not present.
if [ ! -e humans.hdf5 ]
then
    echo "Generating initial conditions for the human mobility box example..."
    # python3 makeIC.py -t gas -f humans.hdf5  # -t particles
    python3 makeIC.py -t gas -n 100 -b 10000 -f humans.hdf5  # -t particles
fi

# Generate acceleration field for river
if [ ! -e river.hdf5 ]
then
    echo "Generating acceleration field for the river..."
    python3 makeRandomRiver.py -b 10000 -g 1000 -f river.hdf5
fi

# Ensure all nodes can access the HDF5 files
#sync
#sleep 2  # Give filesystem time to sync
# ulimit -s unlimited

# Run SWIFT with SLURM environment variables
# SLURM_NTASKS = total number of MPI tasks
# SLURM_CPUS_PER_TASK = number of threads per MPI task
mpirun -np $SLURM_NTASKS swift_mpi --threads=$SLURM_CPUS_PER_TASK \
    -A -s -g -G \
    --hm-river \
    --hm-randomwalk \
    -n 10000 \
    humanMobility.yml

# gdb --args swift -g --threads=4 -n 10000 humanMobility.yml # -A -s
# swift --hm-river --hm-randomwalk --threads=8 -n 1000 humanMobility.yml # -A -s -g -G 
# likwid-perfctr -f -C 0 -g MEMREAD swift -A -s -g -G --hm-river --hm-randomwalk --threads=16 -n 100 humanMobility.yml # -A -s -g -G
# likwid-perfctr -a
