#!/bin/bash

#SBATCH --nodes=1                  # Request 1 node
#SBATCH --ntasks=1                # Single task for compilation
#SBATCH --cpus-per-task=16        # Use all cores for parallel make
#SBATCH --mem=64G                 # Memory for compilation
#SBATCH -J sw-build              # Job name
#SBATCH -o sw.out               # Output file
#SBATCH -e sw.err               # Error file
#SBATCH -p cosma5               # COSMA5 partition
#SBATCH -A durham               # Account
#SBATCH -t 0:30:00             # Increased time for compilation

source ~/swift-env/bin/activate

# Adapt the following lines to your system
#LOCAL_LIBRARY_PATH=/home/dmitry/local/lib
#METIS_PATH=/usr/lib/x86_64-linux-gnu
#PARMETIS_PATH=/usr/lib/x86_64-linux-gnu

module purge
# module load gnu_comp/14.1.0
module load intel_comp/2024.2.0 compiler-rt tbb compiler
module load openmpi/5.0.3
module load fftw/3.3.10
module load gsl
module load parmetis/4.0.3-64bit
module load parallel_hdf5/1.14.4
module load sundials/5.8.0_c8_single
module load ffmpeg
module load python

# Set number of threads for make
export MAKEFLAGS="-j$SLURM_CPUS_PER_TASK"

# The main installation script for SWIFT_ABM

echo "############################"
echo "# Running 'make clean' ... #"
echo "############################"
make clean

echo "#############################"
echo "# Running './configure' ... #"
echo "#############################"
#export CFLAGS="-fsanitize=address -g"
./configure --prefix=/cosma/home/do009/dc-niko3/.local/ \
	    LDFLAGS=\
	    --enable-mpi \
		--enable-ipo \
	    --enable-parallel-hdf5 \
		--with-tbbmalloc \
	    --with-parmetis \
	    --with-hydro=abm \
	    --with-hydro-dimension=2 \
		--with-abm=human-mobility \
		--with-ext-potential=human-mobility \
	    --with-hm=all

	    # --with-hm=random-walk


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
