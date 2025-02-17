#!/bin/bash

# Adapt the following lines to your system
LOCAL_LIBRARY_PATH=/home/dmitry/local/lib
METIS_PATH=/usr/lib/x86_64-linux-gnu
PARMETIS_PATH=/usr/lib/x86_64-linux-gnu

# The main installation script for SWIFT_ABM

echo "############################"
echo "# Running 'make clean' ... #"
echo "############################"
make clean

echo "#############################"
echo "# Running './configure' ... #"
echo "#############################"
./configure CC=gcc \
	    LDFLAGS=-L${LOCAL_LIBRARY_PATH} \
	    --enable-mpi \
	    --enable-parallel-hdf5 \
	    --with-metis=${METIS_PATH} \
	    --with-parmetis=${PARMETIS_PATH} \
	    --with-hydro=abm \
	    --with-hydro-dimension=2 \
		--with-abm=human-mobility \
		--with-ext-potential=human-mobility \
	    --with-hm=river

echo "########################################"
echo "# Running 'make -j\$(nproc --all)' ... #"
echo "########################################"
make -j$(nproc --all)

echo "#####################################"
echo "# Running 'ranlib' on libraries ... #"
echo "#####################################"
ranlib src/.libs/libswiftsim.a argparse/.libs/libargparse.a src/.libs/libswiftsim_mpi.a

echo "############################"
echo "# Running 'make' again ... #"
echo "############################"
make

echo "###################################"
echo "# Running 'sudo make install' ... #"
echo "###################################"
sudo make install
