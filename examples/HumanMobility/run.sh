#!/bin/bash

# Default parameters
NUM_HUMANS=1000

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --num-humans)
            NUM_HUMANS="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --num-humans NUM_HUMANS Number of humans per dimension (default: 1000)"
            echo "  -h, --help              Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

echo "Run Configuration:"
echo "  Number of humans: ${NUM_HUMANS}x${NUM_HUMANS}"

# Determine partition suffix based on SLURM partition
partition_suffix=""
if [[ "${SLURM_JOB_PARTITION}" == "cosma" ]]; then
    partition_suffix="-cosma"
fi

# Basenames for this run (without extensions)
HUMANS=humans-rivers-${NUM_HUMANS}${partition_suffix}
RIVERS=river-rivers-${NUM_HUMANS}${partition_suffix}
HUMANMOBILITY=humanMobility

# Automatically generate directory names based on SLURM parameters
nodes=${SLURM_JOB_NUM_NODES:-1}
ntasks=${SLURM_NTASKS:-1}
cpus_per_task=${SLURM_CPUS_PER_TASK:-16}

# Generate suffix based on resources and partition
DATA=data-rivers-${nodes}-${ntasks}-${cpus_per_task}${partition_suffix}
IMAGES=images-rivers-${nodes}-${ntasks}-${cpus_per_task}${partition_suffix}

# Create the data directory if it doesn't exist
if [ -d "${DATA}" ]; then
    rm -rf ${DATA}/*
fi
mkdir -p ${DATA}

# Check if required input files exist
if [ ! -f "${HUMANS}.hdf5" ]; then
    echo "ERROR: Input file ${HUMANS}.hdf5 not found!"
    echo "Please generate input files first with:"
    echo "  ./submit.sh --gen -- -n ${NUM_HUMANS}"
    exit 1
fi

if [ ! -f "${RIVERS}.hdf5" ]; then
    echo "ERROR: Input file ${RIVERS}.hdf5 not found!"
    echo "Please generate input files first with:"
    echo "  ./submit.sh --gen -- -n ${NUM_HUMANS}"
    exit 1
fi

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
        -A -s -g -G --hm-river --hm-randomwalk -n ${NUM_HUMANS} ${HUMANMOBILITY}.yml
else
    echo "Running PARALLEL version ($ntasks MPI ranks)"
    export OMP_NUM_THREADS=$cpus_per_task
    
    mpirun -n $ntasks \
        ${SWIFT_MPI} --threads=$cpus_per_task \
        --task-dumps=10 -v 1 \
        -A -s -g -G --hm-river --hm-randomwalk -n ${NUM_HUMANS} ${HUMANMOBILITY}.yml
fi

echo "Simulation complete. Check output files and SWIFT task logs for analysis in ${DATA}/"