import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import os

plt.style.use("../../tools/stylesheets/mnras.mplstyle")

import multiprocessing
from functools import partial

def process_file(i, min_x, max_x, min_y, max_y):
    filename_hdf5 = "data/humanMobility_%04d.hdf5" % i  # Adjust the filename pattern as needed
    filename_png = "images/humanMobility_%04d.png" % i

    # Check if the PNG file already exists to avoid reprocessing
    if os.path.exists(filename_png):
        print(f"File {filename_png} already exists. Skipping.")
        return

    # Read data from the HDF5 file
    with h5py.File(filename_hdf5, 'r') as sim:

        # Read the simulation data
        # sim = h5py.File("humanMobility_%04d.hdf5" % snap, "r")
        boxSize = sim["/Header"].attrs["BoxSize"][0]
        time = sim["/Header"].attrs["Time"][0]
        # scheme = sim["/HydroScheme"].attrs["Scheme"]
        # kernel = sim["/HydroScheme"].attrs["Kernel function"]
        # neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
        # eta = sim["/HydroScheme"].attrs["Kernel eta"]
        point_mass = sim["/Parameters"].attrs["PointMassPotential:mass"]
        git = sim["Code"].attrs["Git Revision"]

        pos = sim["/PartType1/Coordinates"][:, :]
        x = pos[:, 0] - boxSize / 2
        y = pos[:, 1] - boxSize / 2
        vel = sim["/PartType1/Velocities"][:, :]
        v_norm = np.sqrt(vel[:, 0] ** 2 + vel[:, 1] ** 2)
        # rho = sim["/PartType0/Densities"][:]
        # u = sim["/PartType0/InternalEnergies"][:]
        # S = sim["/PartType0/Entropies"][:]
        # P = sim["/PartType0/Pressures"][:]

        X = pos[:, 0]
        Y = pos[:, 1]
        U = vel[:, 0]
        V = vel[:, 1]

        # Plot the interesting quantities
        plt.figure(figsize=(7, 7 / 1.6))

        # Calculate the magnitude of the velocity vectors
        # vel_mag = np.sqrt(vel[:, 0]**2 + vel[:, 1]**2)
        M = np.hypot(V, U)
        print("Processing %04d, min/max X=(%4.2f, %4.2f), min/max Y (%4.2f, %4.2f), min/max U=(%4.2f, %4.2f), min/max V=(%4.2f, %4.2f)" % \
              (i, np.min(X), np.max(X), np.min(Y), np.max(Y), np.min(U), np.max(U), np.min(V), np.max(V)))
        print("boxsize:", boxSize)

        # Determine a suitable scale for the arrow lengths
        max_M = np.max(M)
        desired_max_arrow_length = 0.02  # Adjust this value as needed
        scale = max_M / desired_max_arrow_length

        # Plotting the velocity vectors using quiver
        # Create the plot
        fig, ax = plt.subplots()
        ax.set_title("Velocity map of human mobility %04d (scale=%07.2d)\n point-mass=%s" % (i, scale, point_mass))

        # Create the quiver plot
        Q = ax.quiver(
            X,              # X positions
            Y,              # Y positions
            U,              # X components of velocity
            V,              # Y components of velocity
            M,              # Color by velocity magnitude
            scale=scale,    # Scale for arrow lengths
            pivot='tip',
            # angles='xy',
            # scale_units='xy',
            # cmap="PuBu",
            # width=0.025,
            # headwidth=30,
            # headlength=50,
            # headaxislength=45,
            # minshaft=20,
            # minlength=1.0
        )

        ax.scatter(X, Y, color='0.5', s=0.2)

        # Add a colorbar to show the magnitude of velocities
        # plt.colorbar(label='Velocity Magnitude')

        # Add a colorbar associated with the quiver plot
        cbar = fig.colorbar(Q, ax=ax, label='Velocity Magnitude')

        qk = ax.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='figure')
        # Add labels and formatting
        # plt.text(
        #     0.97, 0.97, "${\\rm{Velocity~vectors}}$",
        #     ha="right", va="top", transform=plt.gca().transAxes, backgroundcolor="w"
        # )
        plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
        plt.ylabel("${\\rm{Position}}~y$", labelpad=0)
        plt.xlim(min_x, max_x)
        plt.ylim(min_y, max_y)
        # plt.tight_layout()

        plt.savefig(filename_png, dpi=200)
        plt.close()


if __name__ == '__main__':
    
    # Constants for the arguments
    num_files = int(sys.argv[1])   # Adjust to the number of your HDF5 files
    min_x = int(sys.argv[2])
    max_x = int(sys.argv[3])
    min_y = int(sys.argv[4])
    max_y = int(sys.argv[5])

    # Create a partial function with fixed additional arguments
    partial_process_file = partial(process_file, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)

    file_indices = list(range(num_files))

    # Determine the number of CPUs to use
    num_processes = multiprocessing.cpu_count()

    # Create a pool of worker processes
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Map the partial function to the list of file indices
        pool.map(partial_process_file, file_indices)