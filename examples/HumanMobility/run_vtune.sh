#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-vtune
RIVERS=river-vtune
HUMANMOBILITY=humanMobility
DATA=data-vtune
IMAGES=images-vtune

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVERS
envsubst < humanMobility_template.yml > ${HUMANMOBILITY}.yml

SWIFT=/cosma5/data/durham/dc-niko3/.local/bin/swift_intel2025
SWIFT_MPI=/cosma5/data/durham/dc-niko3/.local/bin/swift_mpi_intel2025

# Use SLURM environment variables for configuration
ntasks=${SLURM_NTASKS:-1}
cpus_per_task=${SLURM_CPUS_PER_TASK:-16}

echo "VTune Analysis Configuration:"
echo "  MPI Tasks: $ntasks"
echo "  OpenMP Threads per Task: $cpus_per_task"

# Enable SWIFT's built-in logging
export SWIFT_TASK_DUMPS=1
export SWIFT_MPIUSE_REPORTS=1
export SWIFT_MEMUSE_REPORTS=1

# Run VTune with automatic serial/parallel detection
if [ $ntasks -eq 1 ]; then
    echo "Running VTune on SERIAL SWIFT"
    
    echo "Running VTune hotspots analysis..."
    vtune -collect hotspots -result-dir vtune_hotspots_intranode \
        ${SWIFT} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

    echo "Generating VTune hotspots summary report..."
    vtune -report summary -result-dir vtune_hotspots_intranode > vtune_hotspots_summary.txt

    echo "Running VTune threading analysis..."
    vtune -collect threading -result-dir vtune_threading_intranode \
        ${SWIFT} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

    echo "Generating VTune threading report..."
    vtune -report summary -result-dir vtune_threading_intranode > vtune_threading_summary.txt
else
    echo "Running VTune on PARALLEL SWIFT"
    
    echo "Running VTune MPI hotspots analysis..."
    srun -n $ntasks -c $cpus_per_task \
        vtune -collect hotspots -result-dir vtune_mpi_hotspots_intranode \
        ${SWIFT_MPI} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

    echo "Generating VTune MPI summary report..."
    vtune -report summary -result-dir vtune_mpi_hotspots_intranode > vtune_mpi_hotspots_summary.txt
fi

echo "VTune analysis complete. Results in vtune_*_summary.txt and vtune_*_intranode/ directories"
echo "To view interactive results, use: vtune-gui vtune_*_intranode"