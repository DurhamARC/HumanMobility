#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-inspector
RIVERS=river-inspector
HUMANMOBILITY=humanMobility
DATA=data-inspector
IMAGES=images-inspector

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVERS
envsubst < humanMobility_template.yml > ${HUMANMOBILITY}.yml

SWIFT=/cosma5/data/durham/dc-niko3/.local/bin/swift_intel2025
SWIFT_MPI=/cosma5/data/durham/dc-niko3/.local/bin/swift_mpi_intel2025

# Use SLURM environment variables for configuration
ntasks=${SLURM_NTASKS:-1}
cpus_per_task=${SLURM_CPUS_PER_TASK:-16}

echo "Intel Inspector Analysis Configuration:"
echo "  MPI Tasks: $ntasks"
echo "  OpenMP Threads per Task: $cpus_per_task"

# Enable SWIFT's built-in logging
export SWIFT_TASK_DUMPS=1
export SWIFT_MPIUSE_REPORTS=1
export SWIFT_MEMUSE_REPORTS=1

# Run Inspector with automatic serial/parallel detection
if [ $ntasks -eq 1 ]; then
    echo "Running Intel Inspector on SERIAL SWIFT"
    
    echo "Running Intel Inspector memory error analysis..."
    inspxe-cl -collect mi2 -result-dir inspector_memory_intranode \
        ${SWIFT} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

    echo "Generating Inspector memory report..."
    inspxe-cl -report summary -result-dir inspector_memory_intranode > inspector_memory_summary.txt

    echo "Running Intel Inspector threading error analysis..."
    inspxe-cl -collect ti2 -result-dir inspector_threading_intranode \
        ${SWIFT} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

    echo "Generating Inspector threading report..."
    inspxe-cl -report summary -result-dir inspector_threading_intranode > inspector_threading_summary.txt
else
    echo "Running Intel Inspector on PARALLEL SWIFT"
    
    echo "Running Intel Inspector MPI memory error analysis..."
    srun -n $ntasks -c $cpus_per_task \
        inspxe-cl -collect mi2 -result-dir inspector_mpi_memory_intranode \
        ${SWIFT_MPI} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

    echo "Generating Inspector MPI memory report..."
    inspxe-cl -report summary -result-dir inspector_mpi_memory_intranode > inspector_mpi_memory_summary.txt
fi

echo "Intel Inspector analysis complete. Results in inspector_*_summary.txt and inspector_*_intranode/ directories"
echo "To view interactive results, use: inspxe-gui inspector_*_intranode"