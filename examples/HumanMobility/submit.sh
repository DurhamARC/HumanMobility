#!/bin/bash

# Submit a job to run the Human Mobility example with specified mode and resources
# ./submit.sh --run -N 1 -n 2 -c 64          # Run on COSMA5 (default, you can append "-p cosma5" optionally)
# ./submit.sh --run -N 1 -n 2 -c 64 -p cosma # Run on COSMA (legacy hardware)
# ./submit.sh --gen -b 50000 -g 5000 -n 500  # Generate with custom parameters
# ./submit.sh --vis -d 1-4-16-cosma          # Visualize specific run results
# ./submit.sh --map -d 1-4-16-cosma          # Generate map with specific suffix
# Note: The script can be extended to other partitions like Cosma7/Cosma8

# Default values
nodes=1
ntasks=1
cpus_per_task=64  # Default for COSMA5, will be adjusted for COSMA
partition="cosma5"  # Default to cosma5, currently implemented for cosma5|cosma
mode=""

# Parse command line arguments
mode_args=()  # Array to store mode-specific arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gen|--vis|--map|--run|--likwid|--maqao|--aps|--vtune|--advisor|--inspector|--scorep)
            mode=$1
            shift
            ;;
        # SLURM-specific options only
        -N|--nodes)
            nodes="$2"
            shift 2
            ;;
        -n|--ntasks)
            ntasks="$2"
            shift 2
            ;;
        -c|--cpus-per-task)
            cpus_per_task="$2"
            shift 2
            ;;
        -p|--partition)
            partition="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 --TYPE [SLURM_OPTIONS] [-- MODE_SPECIFIC_OPTIONS]"
            echo ""
            echo "Types:"
            echo "  gen, vis, map, run, likwid, maqao, aps, vtune, advisor, inspector, scorep"
            echo ""
            echo "SLURM Options:"
            echo "  -N, --nodes NODES           Number of nodes (default: 1)"
            echo "  -n, --ntasks NTASKS         Number of MPI tasks (default: 1)"
            echo "  -c, --cpus-per-task CPUS    CPUs per task (default: 64 for cosma5, 16 for cosma)"
            echo "  -p, --partition PARTITION   Partition: cosma5 (default), cosma"
            echo ""
            echo "Mode-specific options should be passed after -- or directly to the mode scripts"
            echo ""
            echo "Examples:"
            echo "  $0 --run -N 2 -n 4 -c 32"
            echo "  $0 --gen -- -b 50000 -g 5000 -n 500"
            echo "  $0 --vis -- -d 1-4-16-cosma"
            exit 0
            ;;
        --)
            # Everything after -- goes to mode-specific script
            shift
            mode_args=("$@")
            break
            ;;
        *)
            # Unknown options go to mode-specific script
            mode_args+=("$1")
            shift
            ;;
    esac
done

if [ -z "$mode" ]; then
    echo "Usage: $0 --TYPE [SLURM_OPTIONS] [-- MODE_SPECIFIC_OPTIONS]"
    echo "Types: gen, vis, map, run, likwid, maqao, aps, vtune, advisor, inspector, scorep"
    exit 1
fi

# Validate partition
case $partition in
    cosma5|cosma)
        ;;
    *)
        echo "Error: Invalid partition '$partition'. Use 'cosma5' or 'cosma'"
        exit 1
        ;;
esac

# Adjust default cpus_per_task based on partition if not explicitly set
if [[ "$partition" == "cosma" ]] && [[ "$cpus_per_task" == "64" ]]; then
    cpus_per_task=16
fi

name=${mode#--}  # Remove leading -- from mode

# Calculate memory per node (estimate: 4GB per core)
mem_per_node=$((cpus_per_task * ntasks / nodes * 4))G

echo "Submitting to partition: $partition"
echo "Using unified job script: job.sh"
echo "Resources: $nodes nodes, $ntasks tasks, $cpus_per_task cores/task"

# Show mode-specific arguments if any
if [ ${#mode_args[@]} -gt 0 ]; then
    echo "Mode-specific arguments: ${mode_args[*]}"
fi

# Submit job with specified resources and partition
sbatch --job-name="hm-$name" \
       --output="hm-$name-%N-rank%t-$nodes-$ntasks-$cpus_per_task-$partition.out" \
       --error="hm-$name-%N-rank%t-$nodes-$ntasks-$cpus_per_task-$partition.err" \
       --partition="$partition" \
       --nodes=$nodes \
       --ntasks=$ntasks \
       --ntasks-per-node=$((ntasks / nodes)) \
       --cpus-per-task=$cpus_per_task \
       --mem=$mem_per_node \
       job.sh $mode "${mode_args[@]}"