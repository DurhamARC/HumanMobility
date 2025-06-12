#!/bin/bash

# Submit a job to run the Human Mobility example with specified mode and resources
# ./submit.sh --run -N 1 -n 2 -c 64          # Run on COSMA5 (default, you can append "-p cosma5" optionally)
# ./submit.sh --run -N 1 -n 2 -c 64 -p cosma # Run on COSMA (legacy hardware)
# Note: The script can be extended to other partitions like Cosma7/Cosma8

# Default values
nodes=1
ntasks=1
cpus_per_task=64
partition="cosma5"  # Default to cosma5, currently implemented for cosma5|cosma
mode=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gen|--vis|--map|--run|--likwid|--maqao|--aps|--vtune|--advisor|--inspector|--scorep)
            mode=$1
            shift
            ;;
        -N)
            nodes="$2"
            shift 2
            ;;
        -n)
            ntasks="$2"
            shift 2
            ;;
        -c)
            cpus_per_task="$2"
            shift 2
            ;;
        -p|--partition)
            partition="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 --TYPE [-N nodes] [-n ntasks] [-c cpus-per-task] [-p partition]"
            echo "Types: gen, vis, map, run, likwid, maqao, aps, vtune, advisor, inspector, scorep"
            echo "Partitions: cosma5 (default), cosma"
            exit 1
            ;;
    esac
done

if [ -z "$mode" ]; then
    echo "Usage: $0 --TYPE [-N nodes] [-n ntasks] [-c cpus-per-task] [-p partition]"
    echo "Types: gen, vis, map, run, likwid, maqao, aps, vtune, advisor, inspector, scorep"
    echo "Partitions: cosma5 (default), cosma"
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

name=${mode#--}  # Remove leading -- from mode

# Calculate memory per node (estimate: 4GB per core)
mem_per_node=$((cpus_per_task * ntasks / nodes * 4))G

echo "Submitting to partition: $partition"
echo "Using unified job script: job.sh"
echo "Resources: $nodes nodes, $ntasks tasks, $cpus_per_task cores/task"

# Submit job with specified resources and partition - now always use job.sh
sbatch --job-name="hm-$name" \
       --output="hm-$name-$partition-rank%t.out" \
       --error="hm-$name-$partition-rank%t.err" \
       --partition="$partition" \
       --nodes=$nodes \
       --ntasks=$ntasks \
       --ntasks-per-node=$((ntasks / nodes)) \
       --cpus-per-task=$cpus_per_task \
       --mem=$mem_per_node \
       job.sh $mode