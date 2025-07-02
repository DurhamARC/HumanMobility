#!/bin/bash

# Default parameters
BOX_SIZE=100000
GRID_SIZE=10000
NUM_HUMANS=1000

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--box-size)
            BOX_SIZE="$2"
            shift 2
            ;;
        -g|--grid-size)
            GRID_SIZE="$2"
            shift 2
            ;;
        -n|--num-humans)
            NUM_HUMANS="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  -b, --box-size BOX_SIZE     Box size in meters (default: 100000)"
            echo "  -g, --grid-size GRID_SIZE   Grid size for acceleration field (default: 10000)"
            echo "  -n, --num-humans NUM_HUMANS Number of humans per dimension (default: 1000)"
            echo "  -h, --help                  Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Determine partition suffix based on SLURM partition
partition_suffix=""
if [[ "${SLURM_JOB_PARTITION}" == "cosma" ]]; then
    partition_suffix="-cosma"
fi

# Basenames for this run (without extensions)
HUMANS=humans-rivers-${NUM_HUMANS}${partition_suffix}
RIVERS=river-rivers-${NUM_HUMANS}${partition_suffix}

echo "Generation Configuration:"
echo "  Box size: ${BOX_SIZE} m"
echo "  Grid size: ${GRID_SIZE} cells"
echo "  Number of humans: ${NUM_HUMANS}x${NUM_HUMANS}"
echo "  Partition: ${SLURM_JOB_PARTITION:-cosma5}"
echo "  Output files: ${HUMANS}.hdf5, ${RIVERS}.hdf5"
echo ""

# Remove previously generated data and images to regenerate the new ones
rm -f ${HUMANS}.hdf5
rm -f ${RIVERS}.hdf5

# Generate acceleration field for river
if [ ! -e ${RIVERS}.hdf5 ]
then
    echo "Generating acceleration field for the river..."
    python3 makeRivers.py -b ${BOX_SIZE} -g ${GRID_SIZE} -f ${RIVERS}.hdf5
fi

# Generate the initial conditions if they are not present.
if [ ! -e ${HUMANS}.hdf5 ]
then
    echo "Generating initial conditions for the human mobility box example..."
    python3 makeIC.py -t gas -n ${NUM_HUMANS} -b ${BOX_SIZE} -f ${HUMANS}.hdf5
fi

echo "Generation complete!"
echo "  River file: ${RIVERS}.hdf5"
echo "  Humans file: ${HUMANS}.hdf5"