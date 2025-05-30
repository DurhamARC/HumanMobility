#!/bin/bash
#SBATCH --ntasks=1             # Total number of MPI tasks (cores) (max 16)
#SBATCH --nodes=1              # Number of nodes
#SBATCH --ntasks-per-node=1    # MPI tasks per node (max 16)
#SBATCH --cpus-per-task=16     # CPU cores per MPI rank
#SBATCH --mem=120G             # Memory per node
#SBATCH -p cosma               # COSMA5 partition
#SBATCH -A durham              # Account
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=lcgk69@durham.ac.uk
#SBATCH --job-name=hm-perf-both
#SBATCH --output=hm-perf-both.out
#SBATCH --error=hm-perf-both.err

# Load required modules
module purge
module load cosma
module load intel_comp/2025.0.1
module load umf compiler-rt tbb compiler mpi
module load python
module load fftw/3.3.10
module load gsl
module load parmetis/4.0.3-64bit
module load parallel_hdf5/1.14.4
module load sundials/5.8.0_c8_single
module load likwid/5.4.1

module list

echo "=== Running MEM benchmark ==="
./run-mem.sh

echo "=== Running FLOPS_DP benchmark ==="
./run-flops-dp.sh