#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 [--run|--vis]"
    exit 1
fi

mode=$1
name=${mode#--}  # Remove leading -- from mode

sbatch --job-name="hm-$name" \
       --output="hm-$name.out" \
       --error="hm-$name.err" \
       job.sh $mode