#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-scorep
RIVERS=river-scorep
HUMANMOBILITY=humanMobility
DATA=data-scorep
IMAGES=images-scorep

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVERS
envsubst < humanMobility_template.yml > ${HUMANMOBILITY}.yml

SWIFT=/cosma5/data/durham/dc-niko3/.local/bin/swift_intel2025
SWIFT_MPI=/cosma5/data/durham/dc-niko3/.local/bin/swift_mpi_intel2025

# Use SLURM environment variables for configuration
ntasks=${SLURM_NTASKS:-1}
cpus_per_task=${SLURM_CPUS_PER_TASK:-16}

echo "Score-P Analysis Configuration:"
echo "  MPI Tasks: $ntasks"
echo "  OpenMP Threads per Task: $cpus_per_task"

# Enable SWIFT's built-in logging
export SWIFT_TASK_DUMPS=1
export SWIFT_MPIUSE_REPORTS=1
export SWIFT_MEMUSE_REPORTS=1

echo "Running Score-P analysis..."

# Set Score-P environment variables
export SCOREP_ENABLE_PROFILING=true
export SCOREP_ENABLE_TRACING=false
export SCOREP_PROFILING_MAX_CALLPATH_DEPTH=30
export SCOREP_TOTAL_MEMORY=1G

# Note: For Score-P to work properly, SWIFT should be compiled with Score-P instrumentation
echo "Warning: For full Score-P analysis, SWIFT should be recompiled with Score-P instrumentation"
echo "Running with runtime instrumentation only..."

# Run Score-P with automatic serial/parallel detection
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