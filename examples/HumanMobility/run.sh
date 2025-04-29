#!/bin/bash
#SBATCH --ntasks=8             # or 32; total number of MPI tasks (cores)
#SBATCH --nodes=1                       # or 1; number of nodes
#SBATCH --ntasks-per-node=8 # or 32; MPI tasks per node
#SBATCH --cpus-per-task=16  # CPU cores per MPI rank
#SBATCH --mem=32G                      # or 120G; memory per node
#SBATCH -J hm-run
#SBATCH -o hm.out
#SBATCH -e hm.err
#SBATCH -p cosma5                       # cosma5 partition
#SBATCH -A durham
#SBATCH -t 0:30:00
#SBATCH --mail-type=END
#SBATCH --mail-user=lcgk69@durham.ac.uk

module purge

module load cosma
module load intel_comp/2024.2.0 compiler-rt tbb compiler
#module load compiler-rt tbb compiler mpi
module load openmpi/5.0.3
module load fftw/3.3.10
module load gsl
module load parmetis/4.0.3-64bit
module load parallel_hdf5/1.14.4
module load sundials/5.8.0_c8_single
module load python

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
    python3 makeIC.py -t gas -n 50 -b 20000 -f humans.hdf5  # -t particles
fi

# Generate acceleration field for river
if [ ! -e river.hdf5 ]
then
    echo "Generating acceleration field for the river..."
    python3 makeRandomRiver.py -f river.hdf5
fi

# Ensure all nodes can access the HDF5 files
#sync
#sleep 2  # Give filesystem time to sync

# ulimit -s unlimited

# Run SWIFT
# gdb --args swift -g --threads=4 -n 10000 humanMobility.yml # -A -s
# swift -A -s -g -G --hm-river --hm-randomwalk --threads=8 -n 50000 humanMobility.yml # -A -s -g -G 
mpirun --bind-to none -np 8 swift_mpi --threads=16 -A -s -g -G --hm-river --hm-randomwalk -n 100000 humanMobility.yml
