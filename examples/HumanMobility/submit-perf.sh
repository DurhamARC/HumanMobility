#!/bin/bash

sbatch --job-name="hm_likwid" \
       --output="hm_likwid.out" \
       --error="hm_likwid.err" \
       run_likwid.sh