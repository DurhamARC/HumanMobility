#!/bin/bash

#SBATCH --ntasks=16             # or 32; total number of MPI tasks (cores)
#SBATCH --nodes=1                       # or 1; number of nodes
#SBATCH --ntasks-per-node=16 # or 32; MPI tasks per node
#SBATCH --cpus-per-task=16  # CPU cores per MPI rank
#SBATCH --mem=32G                      # or 120G; memory per node
#SBATCH -J hm-vis
#SBATCH -o hm.out
#SBATCH -e hm.err
#SBATCH -p cosma5
#SBATCH -A durham
#SBATCH -t 0:30:00
#SBATCH --mail-type=END
#SBATCH --mail-user=lcgk69@durham.ac.uk

module purge
module load cosma
# module load gnu_comp/14.1.0
module load intel_comp/2024.2.0 compiler-rt tbb compiler
module load openmpi/5.0.3
module load fftw/3.3.10
module load gsl
module load parmetis/4.0.3-64bit
module load parallel_hdf5/1.14.4
module load sundials/5.8.0_c8_single

module load ffmpeg

# Increase file descriptor limit
ulimit -n 4096  # Adjust this number as needed

rm images/*

num_files=1001
min_x=0
max_x=20000
min_y=0
max_y=20000
type="gas"  # or "particles"

# Plot the result
python3 plot_velocity_parallel.py ${num_files} ${min_x} ${max_x} ${min_y} ${max_y} ${type}
# python3 generate_GIF.py ${num_files}

# This command sets:
# - A framerate of 20 frames per second.
# - The input pattern to 'humanMobility_%04d.png' starting from 'humanMobility_0000.png'.
# - The total number of frames to 75 (covering humanMobility_0000.png to humanMobility_0074.png).
# - The output video format to H.264 with yuv420p pixel format for broader compatibility.

# ffmpeg -framerate 60 \
#        -start_number 0 \
#        -i images/humanMobility_%04d.png \
#        -frames:v ${num_files} \
#        -y \
#        -c:v libx264 \
#        -pix_fmt yuv420p \
#        video.mp4
