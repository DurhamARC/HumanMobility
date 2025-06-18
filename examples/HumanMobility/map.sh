#!/bin/bash

# Default parameters
NUM_HUMANS=1000
min_x=0
max_x=100000
min_y=0
max_y=100000

# Parse command line arguments
dir_suffix=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --num-humans)
            NUM_HUMANS="$2"
            shift 2
            ;;
        -d|--dir-suffix)
            dir_suffix="$2"
            shift 2
            ;;
        --min-x)
            min_x="$2"
            shift 2
            ;;
        --max-x)
            max_x="$2"
            shift 2
            ;;
        --min-y)
            min_y="$2"
            shift 2
            ;;
        --max-y)
            max_y="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --num-humans NUM_HUMANS     Number of humans per dimension (default: 1000)"
            echo "  -d, --dir-suffix <suffix>   Suffix part of directory path to read data from"
            echo "  --min-x <value>             Minimum x coordinate for map bounds (default: 0)"
            echo "  --max-x <value>             Maximum x coordinate for map bounds (default: 100000)"
            echo "  --min-y <value>             Minimum y coordinate for map bounds (default: 0)"
            echo "  --max-y <value>             Maximum y coordinate for map bounds (default: 100000)"
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

# If no dir_suffix provided, auto-generate from SLURM parameters
if [ -z "$dir_suffix" ]; then
    nodes=${SLURM_JOB_NUM_NODES:-1}
    ntasks=${SLURM_NTASKS:-1}
    cpus_per_task=${SLURM_CPUS_PER_TASK:-16}
    
    # Determine partition suffix
    partition_suffix=""
    if [[ "${SLURM_JOB_PARTITION}" == "cosma" ]]; then
        partition_suffix="-cosma"
    fi
    
    dir_suffix="${nodes}-${ntasks}-${cpus_per_task}${partition_suffix}"
else
    # Need to determine partition suffix for filename generation
    partition_suffix=""
    if [[ "${SLURM_JOB_PARTITION}" == "cosma" ]] || [[ "$dir_suffix" == *"-cosma" ]]; then
        partition_suffix="-cosma"
    fi
fi

# Basenames for this run (without extensions)
RIVERS=river-rivers-${NUM_HUMANS}${partition_suffix}
IMAGES=images-rivers-${dir_suffix}

echo "Map Generation Configuration:"
echo "  Directory suffix: ${dir_suffix}"
echo "  Images directory: ${IMAGES}"
echo "  Bounds: x=[${min_x}, ${max_x}], y=[${min_y}, ${max_y}]"

# Create images directory if it doesn't exist
mkdir -p ${IMAGES}

# Remove previously generated river images
rm -f ${IMAGES}/rivers.png 2>/dev/null

echo "Generating river map..."
echo "  River file: ${RIVERS}.hdf5"
echo "  Output image: ${IMAGES}/rivers.png"
echo "  Bounds: x[${min_x}, ${max_x}], y[${min_y}, ${max_y}]"

# Plot only the rivers
python3 plot_rivers.py ${RIVERS}.hdf5 ${IMAGES}/rivers.png ${min_x} ${max_x} ${min_y} ${max_y}

echo "Map generation complete! Image saved to: ${IMAGES}/rivers.png"