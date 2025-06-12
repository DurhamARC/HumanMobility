#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-rivers-3
RIVERS=river-rivers-3
HUMANMOBILITY=humanMobility

# Automatically generate directory names based on SLURM parameters
nodes=${SLURM_JOB_NUM_NODES:-1}
ntasks=${SLURM_NTASKS:-1}
cpus_per_task=${SLURM_CPUS_PER_TASK:-16}

# Determine partition suffix based on SLURM partition
partition_suffix=""
if [[ "${SLURM_JOB_PARTITION}" == "cosma" ]]; then
    partition_suffix="-cosma"
fi

# Generate suffix based on resources and partition
DATA=data-rivers-${nodes}-${ntasks}-${cpus_per_task}${partition_suffix}
IMAGES=images-rivers-${nodes}-${ntasks}-${cpus_per_task}${partition_suffix}

# Create the data directory if it doesn't exist
mkdir -p ${DATA}

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVERS
envsubst < humanMobility_template.yml > ${DATA}/${HUMANMOBILITY}.yml

# Select SWIFT binaries based on partition
if [[ "${SLURM_JOB_PARTITION}" == "cosma" ]]; then
    SWIFT=/cosma5/data/durham/dc-niko3/.local/bin/swift_cosma
    SWIFT_MPI=/cosma5/data/durham/dc-niko3/.local/bin/swift_mpi_cosma
else
    SWIFT=/cosma5/data/durham/dc-niko3/.local/bin/swift_intel2025
    SWIFT_MPI=/cosma5/data/durham/dc-niko3/.local/bin/swift_mpi_intel2025
fi

# Enable SWIFT's built-in logging
export SWIFT_TASK_DUMPS=1
export SWIFT_MPIUSE_REPORTS=1
export SWIFT_MEMUSE_REPORTS=1

echo "Auto-detecting execution mode:"
echo "  Partition: ${SLURM_JOB_PARTITION:-cosma5}"
echo "  SLURM_JOB_NUM_NODES: $nodes"
echo "  SLURM_NTASKS: $ntasks"
echo "  SLURM_CPUS_PER_TASK: $cpus_per_task"
echo "  Output directory: ${DATA}"
echo "  SWIFT binary: $(basename ${SWIFT})"

# Change to the data directory so all output files are written there
cd ${DATA}

if [ $ntasks -eq 1 ]; then
    echo "Running SERIAL version (1 MPI rank)"
    export OMP_NUM_THREADS=$cpus_per_task
    
    ${SWIFT} --threads=$cpus_per_task \
        --task-dumps=10 -v 1 \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml
else
    echo "Running PARALLEL version ($ntasks MPI ranks)"
    export OMP_NUM_THREADS=$cpus_per_task
    
    mpirun -n $ntasks \
        ${SWIFT_MPI} --threads=$cpus_per_task \
        --task-dumps=10 -v 1 \
        -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml
fi

echo "Simulation complete. Check output files and SWIFT task logs for analysis in ${DATA}/"