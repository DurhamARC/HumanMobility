#!/bin/bash

# Basenames for this run (without extensions)
HUMANS=humans-rivers-3
RIVERS=river-rivers-3
DATA=data-rivers-3
IMAGES=images-rivers-3

# Remove previously generated data and images to regenerate the new ones
rm -f ${HUMANS}.hdf5
rm -f ${RIVERS}.hdf5
rm -f ${DATA}/*
rm -f ${IMAGES}/*

# Generate acceleration field for river
if [ ! -e ${RIVERS}.hdf5 ]
then
    echo "Generating acceleration field for the river..."
    python3 makeRivers.py -b 100000 -g 10000 -f ${RIVERS}.hdf5
fi

# Generate the initial conditions if they are not present.
if [ ! -e ${HUMANS}.hdf5 ]
then
    echo "Generating initial conditions for the human mobility box example..."
    python3 makeIC.py -t gas -n 1000 -b 100000 -f ${HUMANS}.hdf5  # -t particles
fi