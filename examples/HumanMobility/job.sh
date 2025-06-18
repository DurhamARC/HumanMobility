#!/bin/bash
#SBATCH -A durham               # Account
#SBATCH -t 12:00:00             # Time limit hrs:min:sec
#SBATCH --mail-type=END
##SBATCH --mail-user=<INSERT YOUR EMAIL ADDRESS>
# Note: SLURM partition and other parameters are passed from submit.sh

# Check if an argument was provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 [--gen|--vis|--map|--run|--likwid|--maqao|--aps|--vtune|--advisor|--inspector|--scorep] [mode_specific_params]"
    exit 1
fi

mode=$1
shift  # Remove mode from arguments

case $mode in
    --gen|--vis|--map|--run|--likwid|--maqao|--aps|--vtune|--advisor|--inspector|--scorep)
        ;;
    *)
        echo "Invalid argument. Use --gen, --vis, --map, --run, --likwid, --maqao, --aps, --vtune, --advisor, --inspector, or --scorep"
        exit 1
        ;;
esac

# Store remaining arguments as mode-specific parameters
mode_params="$*"

# Detect partition and set hardware description
if [[ "${SLURM_JOB_PARTITION}" == "cosma" ]]; then
    hardware_desc="COSMA (legacy hardware)"
    intel_comp_version="intel_comp/2024.2.0"  # Use older Intel compiler for compatibility
else
    hardware_desc="COSMA5 (modern hardware)"
    intel_comp_version="intel_comp/2025.0.1"  # Use latest Intel compiler
fi

# Display current job configuration
echo "Job Configuration:"
echo "  Mode: $mode"
echo "  Partition: $hardware_desc"
echo "  Nodes: ${SLURM_JOB_NUM_NODES:-1}"
echo "  MPI Tasks: ${SLURM_NTASKS:-1}"
echo "  CPUs per Task: ${SLURM_CPUS_PER_TASK:-16}"
echo "  Memory per Node: ${SLURM_MEM_PER_NODE:-64G}"

# Show mode-specific parameters if provided
if [ -n "$mode_params" ]; then
    case $mode in
        --gen)
            echo "  Generation Parameters: $mode_params"
            ;;
        --vis)
            echo "  Visualization Parameters: $mode_params"
            ;;
        --map)
            echo "  Map Parameters: $mode_params"
            ;;
        --run)
            echo "  Run Parameters: $mode_params"
            ;;
        --likwid|--maqao|--aps|--vtune|--advisor|--inspector|--scorep)
            echo "  Performance Analysis Parameters: $mode_params"
            ;;
    esac
fi

# Common module loading
module purge
module load cosma
module load $intel_comp_version
if [[ "$intel_comp_version" == "intel_comp/2025.0.1" ]]; then
    module load umf
fi
module load compiler-rt tbb compiler mpi
# module load openmpi/5.0.3/

# Set MPI environment variables
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export I_MPI_FABRICS=shm:ofi

module load python
module load fftw/3.3.10
module load gsl
module load parmetis/4.0.3-64bit
module load parallel_hdf5/1.14.4
module load sundials/5.8.0_c8_single

# Load analysis-specific modules and execute the appropriate script
case $mode in
    --gen)
        echo "Running generation with parameters: $mode_params"
        ./gen.sh $mode_params
        ;;
    --vis)
        module load ffmpeg
        echo "Running visualization with parameters: $mode_params"
        ./vis.sh $mode_params
        ;;
    --map)
        echo "Running map generation with parameters: $mode_params"
        ./map.sh $mode_params
        ;;
    --run)
        echo "Running simulation with parameters: $mode_params"
        ./run.sh $mode_params
        ;;
    --likwid)
        # Load likwid profiler
        module load likwid/5.4.1
        ./run_perf.sh --likwid
        ;;
    --maqao)
        # Load Maqao profiler
        module load maqao
        ./run_perf.sh --maqao
        ;;
    --aps)
        if [[ "${SLURM_JOB_PARTITION}" == "cosma" ]]; then
            echo "Warning: Intel APS 2025 may not be compatible with legacy COSMA hardware"
        fi
        # Load APS profiler
        source /cosma/local/intel/oneAPI_2025.0.1/vtune/2025.0/vtune-vars.sh
        ./run_perf.sh --aps
        ;;
    --vtune)
        if [[ "${SLURM_JOB_PARTITION}" == "cosma" ]]; then
            echo "Warning: Intel VTune 2025 may not be compatible with legacy COSMA hardware"
        fi
        # Load VTune profiler
        source /cosma/local/intel/oneAPI_2025.0.1/vtune/2025.0/vtune-vars.sh
        ./run_perf.sh --vtune
        ;;
    --advisor)
        if [[ "${SLURM_JOB_PARTITION}" == "cosma" ]]; then
            echo "Warning: Intel Advisor 2025 may not be compatible with legacy COSMA hardware"
        fi
        # Load Intel Advisor 2025.0
        source /cosma/local/intel/oneAPI_2025.0.1/advisor/2025.0/advisor-vars.sh
        ./run_perf.sh --advisor
        ;;
    --inspector)
        # Load Intel Inspector
        source /cosma/local/intel/oneAPI_2023.2.0/inspector/2023.2.0/inspxe-vars.sh
        ./run_perf.sh --inspector
        ;;
    --scorep)
        # Load Score-P for detailed instrumentation
        module load scorep/8.4
        ./run_perf.sh --scorep
        ;;
esac
