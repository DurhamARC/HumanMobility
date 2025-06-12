#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-maqao
RIVERS=river-maqao
HUMANMOBILITY=humanMobility
DATA=data-maqao
IMAGES=images-maqao

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVERS
envsubst < humanMobility_template.yml > ${HUMANMOBILITY}.yml

SWIFT=/cosma5/data/durham/dc-niko3/.local/bin/swift_intel2025
SWIFT_MPI=/cosma5/data/durham/dc-niko3/.local/bin/swift_mpi_intel2025

# Use SLURM environment variables for configuration
ntasks=${SLURM_NTASKS:-1}
cpus_per_task=${SLURM_CPUS_PER_TASK:-16}

echo "MAQAO Analysis Configuration:"
echo "  MPI Tasks: $ntasks"
echo "  OpenMP Threads per Task: $cpus_per_task"

# Enable SWIFT's built-in MPI logging for detailed analysis
export SWIFT_TASK_DUMPS=1
export SWIFT_MPIUSE_REPORTS=1
export SWIFT_MEMUSE_REPORTS=1

# Run MAQAO with automatic serial/parallel detection
if [ $ntasks -eq 1 ]; then
    echo "Running MAQAO on SERIAL SWIFT"
    maqao oneview -R1 --output-format=all \
        --output-dir="maqao_serial_$(date +%Y-%m-%d_%H-%M-%S)" -- \
        ${SWIFT} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1
else
    echo "Running MAQAO on PARALLEL SWIFT"
    maqao oneview -R1 --output-format=all \
        --mpi-command="srun -n $ntasks -c $cpus_per_task" \
        --envv_OMP_NUM_THREADS="$cpus_per_task" \
        --envv_SWIFT_TASK_DUMPS="1" \
        --envv_SWIFT_MPIUSE_REPORTS="1" \
        --output-dir="maqao_mpi_$(date +%Y-%m-%d_%H-%M-%S)" -- \
        ${SWIFT_MPI} --threads=$cpus_per_task \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1
fi

echo "MAQAO analysis complete. Check maqao_*/ directory for results."
