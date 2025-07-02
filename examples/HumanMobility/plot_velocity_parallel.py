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

def process_file(i, min_x, max_x, min_y, max_y, _type="gas", river_base="river", data_base="data/humanMobility", image_base="images/humanMobility"):
    river_filename = f"{river_base}.hdf5"
    filename_hdf5 = f"{data_base}_{i:04d}.hdf5"
    filename_png = f"{image_base}_{i:04d}.png"

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(filename_png), exist_ok=True)

    with h5py.File(river_filename, 'r') as sim:
        mass = sim["/Header"].attrs["Mass"]
        box_size = sim["/Header"].attrs["BoxSize"]
        river_group = sim["/RiverGeometry"]
        num_rivers = river_group.attrs.get("num_rivers", 1)
        river_geometries = []
        for idx in range(num_rivers):
            grp = river_group[f"river_{idx}"]
            river_geometries.append({
                "left_bank_x": grp["left_bank_x"][:],
                "left_bank_y": grp["left_bank_y"][:],
                "right_bank_x": grp["right_bank_x"][:],
                "right_bank_y": grp["right_bank_y"][:],
                "center_x": grp["centerline_x"][:],
                "center_y": grp["centerline_y"][:],
            })

    if os.path.exists(filename_png):
        print(f"File {filename_png} already exists. Skipping.")
        return

    # Read data from the HDF5 file
    with h5py.File(filename_hdf5, 'r') as sim:
        boxSize = sim["/Header"].attrs["BoxSize"][0]
        time = sim["/Header"].attrs["Time"][0]
        git = sim["Code"].attrs["Git Revision"]
        part_type = "PartType0" if _type == "gas" else "PartType1"
        pos = sim[f"/{part_type}/Coordinates"][:, :]
        x = pos[:, 0] - boxSize / 2
        y = pos[:, 1] - boxSize / 2
        vel = sim[f"/{part_type}/Velocities"][:, :]
        v_norm = np.sqrt(vel[:, 0] ** 2 + vel[:, 1] ** 2)
        if _type == "gas":
            rho = sim["/PartType0/Densities"][:]
            u = sim["/PartType0/InternalEnergies"][:]
            S = sim["/PartType0/Entropies"][:]
            P = sim["/PartType0/Pressures"][:]

        X = pos[:, 0]
        Y = pos[:, 1]
        U = vel[:, 0]
        V = vel[:, 1]
        plt.figure(figsize=(200, 200 / 1.6))
        M = np.hypot(V, U)
        print("Processing %04d, min/max X=(%4.2f, %4.2f), min/max Y (%4.2f, %4.2f), min/max U=(%4.2f, %4.2f), min/max V=(%4.2f, %4.2f)" % \
              (i, np.min(X), np.max(X), np.min(Y), np.max(Y), np.min(U), np.max(U), np.min(V), np.max(V)))
        print("boxsize:", boxSize)
        max_M = np.max(M)
        desired_max_arrow_length = 0.02
        scale = max_M / desired_max_arrow_length
        fig, ax = plt.subplots()
        ax.set_title("Velocity map of human mobility %04d (scale=%07.2d)\n \"river mass\"=%s" % (i, scale, mass))

        # Plot river banks and centerline
        for river in river_geometries:
            ax.plot(river["left_bank_x"], river["left_bank_y"], '-', color='0.8', linewidth=0.5, label=None)
            ax.plot(river["right_bank_x"], river["right_bank_y"], '-', color='0.8', linewidth=0.5, label=None)
            ax.plot(river["center_x"], river["center_y"], '--', color='0.6', linewidth=0.3, label=None)

        # Create the quiver plot
        Q = ax.quiver(
            X, Y, U, V, M, scale=scale, pivot='tip',
        )
        ax.scatter(X, Y, color='0.5', s=0.2)
        cbar = fig.colorbar(Q, ax=ax, label='Velocity Magnitude')
        qk = ax.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='figure')
        ax.legend(loc='upper right', fontsize='x-small')
        plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
        plt.ylabel("${\\rm{Position}}~y$", labelpad=0)
        plt.xlim(min_x, max_x)
        plt.ylim(min_y, max_y)
        # plt.tight_layout()
        plt.savefig(filename_png, dpi=600)
        plt.close()

if __name__ == '__main__':
    # Usage: python plot_velocity_parallel.py num_files min_x max_x min_y max_y type [river_base] [data_base] [image_base]
    num_files = int(sys.argv[1])
    min_x = int(sys.argv[2])
    max_x = int(sys.argv[3])
    min_y = int(sys.argv[4])
    max_y = int(sys.argv[5])
    _type = sys.argv[6]
    river_base = sys.argv[7] if len(sys.argv) > 7 else "river"
    data_base = sys.argv[8] if len(sys.argv) > 8 else "data/humanMobility"
    image_base = sys.argv[9] if len(sys.argv) > 9 else "images/humanMobility"

    if _type not in ["gas", "particles"]:
        raise ValueError("type must be either 'gas' or 'particles'")

    partial_process_file = partial(
        process_file,
        min_x=min_x,
        max_x=max_x,
        min_y=min_y,
        max_y=max_y,
        _type=_type,
        river_base=river_base,
        data_base=data_base,
        image_base=image_base,
    )

    file_indices = list(range(num_files))

    # Determine the number of CPUs to use
    num_processes = multiprocessing.cpu_count()

    # Create a pool of worker processes
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Map the partial function to the list of file indices
        pool.map(partial_process_file, file_indices)