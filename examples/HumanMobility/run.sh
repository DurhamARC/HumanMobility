
#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e humans.hdf5 ]
then
    echo "Generating initial conditions for the human mobility box example..."
    python3 makeIC.py -f humans.hdf5  #make_hdf5.py
fi

# Run SWIFT
swift --hydro --threads=4 -n 1000 HumanMobility.yml

# Check energy conservation and cooling rate
#python3 plotMobility.py
