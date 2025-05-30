#!/bin/bash

sbatch --job-name="hm-perf-both" \
       --output="hm-perf-both.out" \
       --error="hm-perf-both.err" \
       run-perf.sh