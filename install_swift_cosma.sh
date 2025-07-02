#!/bin/bash

#SBATCH --nodes=1                  # Request 1 node
#SBATCH --ntasks=1                # Single task for compilation
#SBATCH --cpus-per-task=16        # Use all cores for parallel make
#SBATCH --mem=64G                 # Memory for compilation
#SBATCH -J sw-build              # Job name
#SBATCH -o sw_cosma.out               # Output file
#SBATCH -e sw_cosma.err               # Error file
#SBATCH -p cosma               # COSMA5 partition
#SBATCH -A durham               # Account
#SBATCH -t 0:30:00             # Increased time for compilation

source ~/swift-env/bin/activate

# Adapt the following lines to your system
#LOCAL_LIBRARY_PATH=/home/dmitry/local/lib
#METIS_PATH=/usr/lib/x86_64-linux-gnu
#PARMETIS_PATH=/usr/lib/x86_64-linux-gnu

module purge
module load cosma
module load intel_comp/2024.2.0 #intel_comp/2025.0.1
module load umf compiler-rt tbb compiler mpi
# module load openmpi/5.0.3
module load python
module load fftw/3.3.10
module load gsl
module load parmetis/4.0.3-64bit
module load parallel_hdf5/1.14.4
module load sundials/5.8.0_c8_single
module load ffmpeg
# module load likwid/5.4.1

# Set number of threads for make
export MAKEFLAGS="-j$SLURM_CPUS_PER_TASK"
# export AR=/cosma/local/intel/oneAPI_2025.0.1/compiler/2025.0/bin/compiler/llvm-ar
# export LD=/cosma/local/intel/oneAPI_2025.0.1/compiler/2025.0/bin/compiler/llvm-link
# export RANLIB=/cosma/local/intel/oneAPI_2025.0.1/compiler/2025.0/bin/compiler/llvm-ranlib

module list

# The main installation script for SWIFT_ABM

echo "############################"
echo "# Running 'make clean' ... #"
echo "############################"
make clean

echo "#############################"
echo "# Running './configure' ... #"
echo "#############################"
#===============================================================================
# SWIFT Configuration for Human Mobility Analysis
#===============================================================================
#
# This configuration builds SWIFT with the following features:
#
# DEBUGGING & PROFILING:
#   --enable-debug           : Enable debug symbols and consistency checks
#   --enable-task-debugging  : Enable detailed task-level performance analysis
#                             (generates thread_info_MPI-step*.dat files)
#
# PARALLELIZATION:
#   --enable-mpi            : MPI parallel support for multi-node execution
#   --enable-parallel-hdf5  : Parallel HDF5 I/O for efficient data output
#   --with-tbbmalloc        : Intel TBB memory allocator for performance
#   --with-parmetis         : ParMETIS library for domain decomposition
#
# HUMAN MOBILITY MODEL:
#   --with-hydro=abm              : Agent-based modeling scheme (a placeholder for future ABM schemes, separating it from the hydro scheme)
#   --with-hydro-dimension=2      : 2D simulation space
#   --with-abm=human-mobility     : Human mobility implementation (a placeholder for future ABM implementations)
#   --with-ext-potential=human-mobility : External potential for human mobility
#
# HUMAN MOBILITY SCENARIOS (--with-hm options):
#   --with-hm=none          : No human mobility scenario (default)
#                             Disables all HM-specific features
#
#   --with-hm=all           : Enable all human mobility scenarios
#                             Includes both river and random-walk cases
#                             Activates HM_CASE_RIVER and HM_CASE_RANDOMWALK
#
#   --with-hm=river         : River-based human mobility scenario
#                             Enables river geography
#                             Activates HM_CASE_RIVER preprocessor definition
#
#   --with-hm=random-walk   : Random walk human mobility scenario  
#                             Initiate human random walk which consequently is influenced by other humans and geographic features
#                             Activates HM_CASE_RANDOMWALK preprocessor definition
#
# NOTE: Human mobility scenarios require:
#       - --with-abm=human-mobility (mandatory)
#       - --with-hydro=abm (mandatory)
#       - --with-ext-potential=human-mobility (recommended)
#       - --with-hydro-dimension=2 (typical for HM studies)
#
#===============================================================================
#/configure --prefix=/cosma/home/do009/dc-niko3/.local/ \
#export CFLAGS="-fsanitize=address -g"
# --program-suffix=_intel2025 \
./configure --prefix=/cosma5/data/durham/dc-niko3/.local/ \
    --program-suffix=_cosma \
    CFLAGS="-Wno-error=gnu-folding-constant" \
    LDFLAGS= \
    --enable-debug \
    --enable-task-debugging \
    --enable-mpi \
    --enable-parallel-hdf5 \
    --with-tbbmalloc \
    --with-parmetis \
    --with-hydro=abm \
    --with-hydro-dimension=2 \
    --with-abm=human-mobility \
    --with-ext-potential=human-mobility \
    --with-hm=all

echo "########################################"
echo "# Running 'make -j\$(nproc --all)' ... #"
echo "########################################"
make

echo "#####################################"
echo "# Running 'ranlib' on libraries ... #"
echo "#####################################"
ranlib src/.libs/libswiftsim.a argparse/.libs/libargparse.a src/.libs/libswiftsim_mpi.a

echo "############################"
echo "# Running 'make' again ... #"
echo "############################"
make

echo "###################################"
echo "# Running 'sudo make install' ... #"
echo "###################################"
make install
