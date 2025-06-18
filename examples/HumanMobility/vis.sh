#!/bin/bash

# Default parameters
dir_suffix=""
num_files=18
min_x=4000
max_x=6000
min_y=4000
max_y=6000
type="gas"
NUM_HUMANS=1000

# Parse command line arguments
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
        --num-files)
            num_files="$2"
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
        --type)
            type="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --num-humans NUM_HUMANS     Number of humans per dimension (default: 1000)"
            echo "  -d, --dir-suffix <suffix>   Suffix part of directory path to read data from"
            echo "  --num-files <num>           Number of output files to process (default: 18)"
            echo "  --min-x <value>             Minimum x coordinate for visualization (default: 4000)"
            echo "  --max-x <value>             Maximum x coordinate for visualization (default: 6000)"
            echo "  --min-y <value>             Minimum y coordinate for visualization (default: 4000)"
            echo "  --max-y <value>             Maximum y coordinate for visualization (default: 6000)"
            echo "  --type <type>               Visualization type: gas or particles (default: gas)"
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
HUMANS=humans-rivers-${NUM_HUMANS}${partition_suffix}
RIVERS=river-rivers-${NUM_HUMANS}${partition_suffix}
HUMANMOBILITY=humanMobility
DATA=data-rivers-${dir_suffix}
IMAGES=images-rivers-${dir_suffix}

echo "Visualization Configuration:"
echo "  Directory suffix: ${dir_suffix}"
echo "  Data directory: ${DATA}"
echo "  Images directory: ${IMAGES}"
echo "  Number of files: ${num_files}"
echo "  Visualization bounds: x=[${min_x}, ${max_x}], y=[${min_y}, ${max_y}]"
echo "  Visualization type: ${type}"

# Create images directory if it doesn't exist
mkdir -p ${IMAGES}

# Remove existing images
rm -f ${IMAGES}/* # don't remore if locally

# Check if data directory exists
if [ ! -d "${DATA}" ]; then
    echo "Error: Data directory '${DATA}' not found!"
    echo "Make sure you have run a simulation first or specify the correct directory suffix with -d"
    exit 1
fi

echo "Processing ${num_files} output files..."

# Use the new configurable parameters
python3 plot_velocity_parallel.py ${num_files} ${min_x} ${max_x} ${min_y} ${max_y} ${type} ${RIVERS} ${DATA}/${HUMANMOBILITY} ${IMAGES}/${HUMANMOBILITY}

echo "Visualization complete! Images saved to: ${IMAGES}/"

# Optional: Generate video (commented out by default)
# python3 generate_GIF.py ${num_files}

# This command sets:
# - A framerate of 20 frames per second.
# - The input pattern to 'humanMobility_%04d.png' starting from 'humanMobility_0000.png'.
# - The total number of frames to 75 (covering humanMobility_0000.png to humanMobility_0074.png).
# - The output video format to H.264 with yuv420p pixel format for broader compatibility.

# ffmpeg -framerate 100 \ # uncomment this line if locally
#        -start_number 0 \
#        -i ${IMAGES}/humanMobility_%04d.png \
#        -frames:v ${num_files} \
#        -y \
#        -c:v libx264 \
#        -pix_fmt yuv420p \
#        ${IMAGES}/video.mp4
