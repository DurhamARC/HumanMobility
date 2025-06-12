#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-advisor
RIVERS=river-advisor
HUMANMOBILITY=humanMobility
DATA=data-advisor
IMAGES=images-advisor

# Render the YAML config from template
export DATA HUMANMOBILITY HUMANS RIVERS
envsubst < humanMobility_template.yml > ${HUMANMOBILITY}.yml

SWIFT=/cosma5/data/durham/dc-niko3/.local/bin/swift_intel2025
SWIFT_MPI=/cosma5/data/durham/dc-niko3/.local/bin/swift_mpi_intel2025

# Use SLURM environment variables for configuration
ntasks=${SLURM_NTASKS:-1}
cpus_per_task=${SLURM_CPUS_PER_TASK:-16}

echo "Intel Advisor Analysis Configuration:"
echo "  MPI Tasks: $ntasks"
echo "  OpenMP Threads per Task: $cpus_per_task"

# Enable SWIFT's built-in logging
export SWIFT_TASK_DUMPS=1
export SWIFT_MPIUSE_REPORTS=1
export SWIFT_MEMUSE_REPORTS=1

# Run Advisor with automatic serial/parallel detection
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