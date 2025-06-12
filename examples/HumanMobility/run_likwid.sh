#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-likwid
RIVERS=river-likwid
HUMANMOBILITY=humanMobility
DATA=data-likwid
IMAGES=images-likwid

# Create the data directory if it doesn't exist
mkdir -p ${DATA}

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVERS
envsubst < humanMobility_template.yml > ${DATA}/${HUMANMOBILITY}.yml

SWIFT=/cosma5/data/durham/dc-niko3/.local/bin/swift_intel2025
SWIFT_MPI=/cosma5/data/durham/dc-niko3/.local/bin/swift_mpi_intel2025

# Use SLURM environment variables for configuration
ntasks=${SLURM_NTASKS:-1}
cpus_per_task=${SLURM_CPUS_PER_TASK:-16}

echo "LIKWID Analysis Configuration:"
echo "  MPI Tasks: $ntasks"
echo "  OpenMP Threads per Task: $cpus_per_task"
echo "  Output directory: ${DATA}"

# Change to the data directory so all output files are written there
cd ${DATA}

# Enable SWIFT's built-in logging
export SWIFT_TASK_DUMPS=1
export SWIFT_MPIUSE_REPORTS=1
export SWIFT_MEMUSE_REPORTS=1

# Run LIKWID with automatic serial/parallel detection
if [ $ntasks -eq 1 ]; then
    echo "Running LIKWID on SERIAL SWIFT"
    
    echo "=== Running LIKWID MEM benchmark ==="
    likwid-perfctr -C 0-$((cpus_per_task-1)) -g MEM \
        ${SWIFT} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1 \
        > likwid_mem_summary.txt 2>&1

    echo "=== Running LIKWID FLOPS_DP benchmark ==="
    likwid-perfctr -C 0-$((cpus_per_task-1)) -g FLOPS_DP \
        ${SWIFT} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1 \
        > likwid_flops_summary.txt 2>&1

    echo "=== Running LIKWID L3 Cache benchmark ==="
    likwid-perfctr -C 0-$((cpus_per_task-1)) -g L3 \
        ${SWIFT} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1 \
        > likwid_cache_summary.txt 2>&1
else
    echo "Running LIKWID on PARALLEL SWIFT"
    
    echo "=== Running LIKWID MEM benchmark (MPI) ==="
    srun -n $ntasks -c $cpus_per_task \
        likwid-mpirun -np $ntasks -g MEM \
        ${SWIFT_MPI} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1 \
        > likwid_mpi_mem_summary.txt 2>&1

    echo "=== Running LIKWID FLOPS_DP benchmark (MPI) ==="
    srun -n $ntasks -c $cpus_per_task \
        likwid-mpirun -np $ntasks -g FLOPS_DP \
        ${SWIFT_MPI} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1 \
        > likwid_mpi_flops_summary.txt 2>&1
fi

echo "LIKWID analysis complete. Results in likwid_*_summary.txt files"
echo "Check the generated summary files for detailed performance metrics"