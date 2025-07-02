#!/bin/bash

# Unified performance analysis script for SWIFT Human Mobility simulations
# Usage: ./run_perf.sh --TOOL
# where TOOL can be: likwid, maqao, aps, vtune, advisor, inspector, scorep

# Check if an argument was provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 [--likwid|--maqao|--aps|--vtune|--advisor|--inspector|--scorep]"
    exit 1
fi

mode=$1

case $mode in
    --likwid|--maqao|--aps|--vtune|--advisor|--inspector|--scorep)
        ;;
    *)
        echo "Invalid argument. Use --likwid, --maqao, --aps, --vtune, --advisor, --inspector, or --scorep"
        exit 1
        ;;
esac

# Extract tool name from argument
TOOL=${mode#--}

# Basenames for this run (without extensions)
HUMANS=humans-${TOOL}
RIVERS=river-${TOOL}
HUMANMOBILITY=humanMobility

# Automatically generate directory names based on SLURM parameters
nodes=${SLURM_JOB_NUM_NODES:-1}
ntasks=${SLURM_NTASKS:-1}
cpus_per_task=${SLURM_CPUS_PER_TASK:-16}

# Determine partition suffix
partition_suffix=""
if [[ "${SLURM_JOB_PARTITION}" == "cosma" ]]; then
    partition_suffix="-cosma"
fi

# Generate suffix based on resources and tool
DATA=data-${TOOL}-${nodes}-${ntasks}-${cpus_per_task}${partition_suffix}
IMAGES=images-${TOOL}-${nodes}-${ntasks}-${cpus_per_task}${partition_suffix}

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

echo "${TOOL^^} Analysis Configuration:"
echo "  Tool: $TOOL"
echo "  Partition: ${SLURM_JOB_PARTITION:-cosma5}"
echo "  MPI Tasks: $ntasks"
echo "  OpenMP Threads per Task: $cpus_per_task"
echo "  Output directory: ${DATA}"

# Change to the data directory so all output files are written there
cd ${DATA}

# Enable SWIFT's built-in logging (disable for APS to avoid conflicts)
if [[ "$TOOL" != "aps" ]]; then
    export SWIFT_TASK_DUMPS=1
    export SWIFT_MPIUSE_REPORTS=1
    export SWIFT_MEMUSE_REPORTS=1
fi

# Tool-specific analysis functions
run_likwid() {
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
}

run_maqao() {
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
}

run_aps() {
    if [ $ntasks -eq 1 ]; then
        echo "Running Intel APS on SERIAL SWIFT"
        
        echo "Running Intel APS performance snapshot analysis..."
        aps -r aps_serial_intranode \
            ${SWIFT} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml

        echo "Generating APS summary report..."
        aps --report=summary --result-dir=aps_serial_intranode > aps_serial_summary.txt

        echo "Generating APS detailed report..."
        aps --report=detailed --result-dir=aps_serial_intranode > aps_serial_detailed.txt
    else
        echo "Running Intel APS on PARALLEL SWIFT"
        
        echo "Running Intel APS MPI performance snapshot analysis..."
        srun -n $ntasks -c $cpus_per_task \
            aps -r aps_mpi_intranode \
            ${SWIFT_MPI} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml

        echo "Generating APS MPI summary report..."
        aps --report=summary --result-dir=aps_mpi_intranode > aps_mpi_summary.txt

        echo "Generating APS MPI detailed report..."
        aps --report=detailed --result-dir=aps_mpi_intranode > aps_mpi_detailed.txt
    fi
    echo "Intel APS analysis complete. Results in aps_*_summary.txt and aps_*_detailed.txt files"
	echo "Result directories: aps_*_intranode/"
	echo ""
	echo "Key APS metrics to check:"
	echo "  - CPU utilization and efficiency"
	echo "  - Memory bandwidth utilization"
	echo "  - Elapsed time and performance characteristics"
	echo "  - MPI communication overhead (for parallel runs)"
	echo ""
	echo "To view results:"
	echo "  cat aps_*_summary.txt     # Quick overview"
	echo "  cat aps_*_detailed.txt    # Detailed analysis"
}

run_vtune() {
    if [ $ntasks -eq 1 ]; then
        echo "Running VTune on SERIAL SWIFT"
        
        echo "Running VTune hotspots analysis..."
        vtune -collect hotspots -result-dir vtune_hotspots_intranode \
            ${SWIFT} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

        echo "Generating VTune hotspots summary report..."
        vtune -report summary -result-dir vtune_hotspots_intranode > vtune_hotspots_summary.txt

        echo "Running VTune threading analysis..."
        vtune -collect threading -result-dir vtune_threading_intranode \
            ${SWIFT} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

        echo "Generating VTune threading report..."
        vtune -report summary -result-dir vtune_threading_intranode > vtune_threading_summary.txt
    else
        echo "Running VTune on PARALLEL SWIFT"
        
        echo "Running VTune MPI hotspots analysis..."
        srun -n $ntasks -c $cpus_per_task \
            vtune -collect hotspots -result-dir vtune_mpi_hotspots_intranode \
            ${SWIFT_MPI} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

        echo "Generating VTune MPI summary report..."
        vtune -report summary -result-dir vtune_mpi_hotspots_intranode > vtune_mpi_hotspots_summary.txt
    fi
    echo "VTune analysis complete. Results in vtune_*_summary.txt and vtune_*_intranode/ directories"
	echo "To view interactive results, use: vtune-gui vtune_*_intranode"
}

run_advisor() {
    if [ $ntasks -eq 1 ]; then
        echo "Running Intel Advisor on SERIAL SWIFT"
        
        echo "Running Intel Advisor survey analysis..."
        advisor --collect=survey --project-dir=advisor_intranode \
            ${SWIFT} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

        echo "Generating Advisor survey report..."
        advisor --report=survey --project-dir=advisor_intranode --format=text > advisor_survey_report.txt

        echo "Running Intel Advisor tripcounts analysis..."
        advisor --collect=tripcounts --project-dir=advisor_intranode \
            ${SWIFT} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

        echo "Generating Advisor tripcounts report..."
        advisor --report=tripcounts --project-dir=advisor_intranode --format=text > advisor_tripcounts_report.txt

        echo "Running Intel Advisor map analysis for vectorization opportunities..."
        advisor --collect=map --project-dir=advisor_intranode \
            ${SWIFT} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

        echo "Generating Advisor map report..."
        advisor --report=map --project-dir=advisor_intranode --format=text > advisor_map_report.txt
    else
        echo "Running Intel Advisor on PARALLEL SWIFT"
        
        echo "Running Intel Advisor MPI survey analysis..."
        srun -n $ntasks -c $cpus_per_task \
            advisor --collect=survey --project-dir=advisor_mpi_intranode \
            ${SWIFT_MPI} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

        echo "Generating Advisor MPI survey report..."
        advisor --report=survey --project-dir=advisor_mpi_intranode --format=text > advisor_mpi_survey_report.txt
    fi
    echo "Intel Advisor analysis complete. Results in advisor_*_report.txt and advisor_*_intranode/ directory"
	echo "To view interactive results, use: advisor-gui advisor_*_intranode"
}

run_inspector() {
    if [ $ntasks -eq 1 ]; then
        echo "Running Intel Inspector on SERIAL SWIFT"
        
        echo "Running Intel Inspector memory error analysis..."
        inspxe-cl -collect mi2 -result-dir inspector_memory_intranode \
            ${SWIFT} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

        echo "Generating Inspector memory report..."
        inspxe-cl -report summary -result-dir inspector_memory_intranode > inspector_memory_summary.txt

        echo "Running Intel Inspector threading error analysis..."
        inspxe-cl -collect ti2 -result-dir inspector_threading_intranode \
            ${SWIFT} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

        echo "Generating Inspector threading report..."
        inspxe-cl -report summary -result-dir inspector_threading_intranode > inspector_threading_summary.txt
    else
        echo "Running Intel Inspector on PARALLEL SWIFT"
        
        echo "Running Intel Inspector MPI memory error analysis..."
        srun -n $ntasks -c $cpus_per_task \
            inspxe-cl -collect mi2 -result-dir inspector_mpi_memory_intranode \
            ${SWIFT_MPI} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1

        echo "Generating Inspector MPI memory report..."
        inspxe-cl -report summary -result-dir inspector_mpi_memory_intranode > inspector_mpi_memory_summary.txt
    fi
    echo "Intel Inspector analysis complete. Results in inspector_*_summary.txt and inspector_*_intranode/ directories"
	echo "To view interactive results, use: inspxe-gui inspector_*_intranode"
}

run_scorep() {
    echo "Running Score-P analysis..."

    # Set Score-P environment variables
    export SCOREP_ENABLE_PROFILING=true
    export SCOREP_ENABLE_TRACING=false
    export SCOREP_PROFILING_MAX_CALLPATH_DEPTH=30
    export SCOREP_TOTAL_MEMORY=1G

	# Note: For Score-P to work properly, SWIFT should be compiled with Score-P instrumentation
    echo "Warning: For full Score-P analysis, SWIFT should be recompiled with Score-P instrumentation"
    echo "Running with runtime instrumentation only..."

    if [ $ntasks -eq 1 ]; then
        echo "Running Score-P on SERIAL SWIFT"
        
        ${SWIFT} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1
    else
        echo "Running Score-P on PARALLEL SWIFT"
        
        srun -n $ntasks -c $cpus_per_task \
            ${SWIFT_MPI} --threads=$cpus_per_task \
            -A -s -g -G --hm-river --hm-randomwalk -n 1000 ${HUMANMOBILITY}.yml --task-dumps=1
    fi

    # Generate Score-P report if profile data exists
    if ls scorep-* 1> /dev/null 2>&1; then
        echo "Generating Score-P summary report..."
        scorep-score scorep-*/profile.cubex > scorep_intranode_summary.txt
        
        echo "Score-P analysis complete. Results in scorep_intranode_summary.txt and scorep-*/ directory"
        echo "Use 'cube scorep-*/profile.cubex' for interactive analysis"
    else
        echo "No Score-P profile data generated. Consider recompiling SWIFT with Score-P instrumentation."
        echo "To enable full Score-P analysis, recompile SWIFT with:"
        echo "  CC='scorep-gcc' CXX='scorep-g++' ./configure [options]"
    fi
}

# Execute the appropriate analysis function
case $TOOL in
    likwid)
        run_likwid
        ;;
    maqao)
        run_maqao
        ;;
    aps)
        run_aps
        ;;
    vtune)
        run_vtune
        ;;
    advisor)
        run_advisor
        ;;
    inspector)
        run_inspector
        ;;
    scorep)
        run_scorep
        ;;
esac

echo "Check the generated files for detailed performance metrics in ${DATA}/"