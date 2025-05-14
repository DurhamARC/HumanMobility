#!/bin/bash

#SBATCH --ntasks=1             # Total number of MPI tasks (cores) (max 16)
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1    # MPI tasks per node (max 16)
#SBATCH --cpus-per-task=16      # CPU cores per MPI rank
#SBATCH --mem=120G               # Memory per node
#SBATCH -p cosma               # COSMA5 partition
#SBATCH -A durham               # Account
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=lcgk69@durham.ac.uk

# Check if an argument was provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 [--run|--vis]"
    exit 1
fi

mode=$1

case $mode in
    --run|--vis)
        ;;
    *)
        echo "Invalid argument. Use --run or --vis"
        exit 1
        ;;
esac

# Common module loading for both modes
module purge
module load cosma
module load intel_comp/2025.0.1 # gnu_comp/14.1.0 intel_comp/2024.2.0
module load umf compiler-rt tbb compiler mpi
# module load openmpi/5.0.3/
module load python
module load fftw/3.3.10
module load gsl
module load parmetis/4.0.3-64bit
module load parallel_hdf5/1.14.4
module load sundials/5.8.0_c8_single
module load likwid/5.4.1


# Additional module for visualisation
if [ "$mode" = "--vis" ]; then
    module load ffmpeg
fi

module list

# Execute the appropriate script
case $mode in
    --run)
        ./run.sh
        ;;
    --vis)
        ./visualise.sh
        ;;
esac
