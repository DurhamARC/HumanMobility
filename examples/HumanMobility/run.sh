
#!/bin/bash

# Remove previously generated data and images to regenerate the new ones
rm humans.hdf5
rm data/*
rm images/*

# Generate the initial conditions if they are not present.
if [ ! -e humans.hdf5 ]
then
    echo "Generating initial conditions for the human mobility box example..."
    python3 makeIC.py -f humans.hdf5  #make_hdf5.py
fi

# Run SWIFT
swift --hydro --threads=4 -n 1000 humanMobility.yml

# Plot the result
python3 plot_velocity_parallel.py 450 # number of snapshots

python3 generate_GIF.py
